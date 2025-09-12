#!/usr/bin/env python3

"""
Vyvanse (LDX → d-amphetamine) + Dex IR — PK & perceived (PD) overlay

• PK curves (solid): model plasma concentration. Vyvanse and Dex share the same
  elimination half-life (ke), while Dex absorbs faster (higher ka).
• Perceived curves (dotted): perceived/pharmacodynamic response derived from PK
  via a shaping kernel (biexponential band-pass). This helps illustrate why
  perceived effect can peak later and wear off sooner than plasma.
• Stop-after projections: branches that show what happens if no further Dex doses
  are taken after a given point (auto-generated for each Dex dose).

Usage:
  - Adjust the VYVANSE list below for Vyvanse dose(s). Each entry is (time_of_day, dose_mg).
  - Adjust the DEX list below for Dex dose(s). Each entry is (time_of_day, dose_mg).
  - Any number of Dex doses are supported; labels and projections update automatically.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import math
from utils.save_utils import save_figure_safely
from utils.dosing_utils import (
    shots_to_caffeine_mg,
    aeropress_scoops_to_caffeine_mg,
    grams_to_caffeine_mg,
    vyvanse_cap_to_dex_eq,
)
from utils.pk_models import bateman, caffeine_total_curve, curves_from_schedule
from utils.plot_utils import label_hour
from utils.style import DEX_BASE_COLORS, COLORS

 

# === Perceived-effect kernel ===

# === Global defaults ===
RES_MIN = 1  # minutes per sample
KA_VYV, KA_DEX = 0.80, 1.00
KE_AMP = np.log(2) / 11.0
DEFAULT_OUTPUT = "graph-vyvanse-dex-pk-vs-perceived.svg"

# Example dosing schedule 
VYVANSE = [(8.0, vyvanse_cap_to_dex_eq(30.0))]
DEX = [(8.0, 5.0), (11.0, 5.0), (13.0, 5.0)]  # Dex 5mg at 8am, 11am, 1pm

# Caffeine schedule examples (each item: (time_h, caffeine_mg[, drink_duration_min[, parts]]))
# - If duration is provided and > 0, the dose is split uniformly across that duration (in minutes).
# - Optional 4th field `parts` lets you set how many equal sub-doses to split into across that time.
# - Use helper conversions to estimate mg: shots_to_caffeine_mg(), aeropress_scoops_to_caffeine_mg(), grams_to_caffeine_mg().
#
# Example: 2-shot espresso at 9am sipped over 20 minutes
# CAFFEINE = [(9.0, shots_to_caffeine_mg(2), 20)]
#
# Example: split the same 20 minutes into 4 equal parts explicitly
# CAFFEINE = [(9.0, shots_to_caffeine_mg(2), 20, 4)]
#
# Example: AeroPress with 1 scoops at 1pm in 60 minutes split into 4 parts
# CAFFEINE = [(13.0, aeropress_scoops_to_caffeine_mg(1.0), 60, 4)]
CAFFEINE = [(13.0, aeropress_scoops_to_caffeine_mg(1.0), 60)]

# === Time window helpers ===
def compute_time_window(vyv_schedule, dex_schedule, caffeine_schedule=None, default_start=8.0):
    times = [td for td, _ in (vyv_schedule or [])] + [td for td, _ in (dex_schedule or [])]
    if caffeine_schedule:
        for item in caffeine_schedule:
            if isinstance(item, (list, tuple)) and len(item) >= 1:
                times.append(float(item[0]))
    if times:
        first = min(times)
        start = float(math.floor(first))
    else:
        start = float(default_start)
    end = start + 24.0
    return start, end

# Perceived-effect (PD) parameters
default_PD = dict(
    pd_floor=0.05,     # mask PD below this threshold for clarity
    # Kernel parameters
    dex_tau_r=0.5,     # Dex perceived rise time constant (h)
    dex_tau_d=3.0,     # Dex perceived decay time constant (h)
    dex_gain=1.0,      # Dex PD gain (pre-matching)
    vyv_tau_r=1.0,     # Vyvanse perceived rise time constant (h)
    vyv_tau_d=6.0,     # Vyvanse perceived decay time constant (h)
    vyv_gain=1.0,      # Vyvanse PD gain (pre-matching)
    pd_peak_scale=1.0, # Target PD_peak = this * PK_peak
    pd_max_scale=1.1   # Clamp PD ≤ this * PK_peak
)

# === Build PK curves ===
def build_pk(t):
    vyv_curves = curves_from_schedule(t, VYVANSE, KA_VYV, KE_AMP)
    dex_curves = curves_from_schedule(t, DEX, KA_DEX, KE_AMP)
    vyv_sum = np.nansum(np.vstack([np.nan_to_num(c) for c in vyv_curves]), axis=0) if vyv_curves else np.zeros_like(t)
    dex_sum = np.nansum(np.vstack([np.nan_to_num(c) for c in dex_curves]), axis=0) if dex_curves else np.zeros_like(t)
    return vyv_sum, vyv_curves, dex_curves, vyv_sum + dex_sum

# === PK→PD transform ===
def mask_below(arr, thresh):
    """Return a copy with values < thresh masked as NaN for plotting clarity."""
    a = np.array(arr, copy=True)
    m = a < thresh
    a[m] = np.nan
    return a

def pd_kernel_biexp(dt_h, tau_r, tau_d, gain=1.0):
    """Zero-DC biexponential kernel for PD shaping (band-pass like).

    k(t) = e^{-t/tau_r}/tau_r - e^{-t/tau_d}/tau_d, t>=0
    Normalized so the positive lobe area equals 1 (stable amplitude),
    then scaled by `gain`.
    """
    tau_r = max(1e-3, float(tau_r))
    tau_d = max(tau_r + 1e-3, float(tau_d))
    t_max = int(np.ceil(8 * tau_d / dt_h))
    t_k = np.arange(0, t_max + 1) * dt_h
    k = np.exp(-t_k / tau_r) / tau_r - np.exp(-t_k / tau_d) / tau_d
    # Normalize by positive lobe area (L1+), more robust than peak norm
    pos_area = np.sum(np.clip(k, 0, None)) * dt_h
    if pos_area <= 1e-9:
        pos_area = 1.0
    k = (k / pos_area) * gain
    return k

def apply_pd_kernel(curve, dt_h, tau_r, tau_d, gain=1.0, match='peak', peak_scale=1.0, clamp_scale=None):
    """Apply kernel PD and scale to match PK peak proportion.

    - match='peak': rescale PD so PD_peak = peak_scale * PK_peak.
    - clamp_scale: if set (e.g., 1.2), cap PD at clamp_scale * PK_peak.
    """
    c = np.nan_to_num(curve)
    k = pd_kernel_biexp(dt_h, tau_r, tau_d, gain)
    y = np.convolve(c, k, mode='full')[:len(c)]
    y = np.maximum(y, 0.0)
    if match == 'peak':
        pk_peak = float(np.max(c)) if c.size else 1.0
        pd_peak = float(np.max(y)) if y.size else 1.0
        if pd_peak > 1e-9:
            y = y * ((peak_scale * pk_peak) / pd_peak)
        if clamp_scale is not None and pk_peak > 0:
            y = np.minimum(y, clamp_scale * pk_peak)
    return y

# === Plotting ===
def plot_overlay(t, vyv_sum, vyv_pk_curves, dex_pk_curves, total_pk, PD, t_start, t_end):
    fig, ax = plt.subplots(figsize=(13, 7))
    # Colors come from centralized DEX_BASE_COLORS / COLORS

    # Helper to mask pre-dose zero baselines for clarity
    def mask_before(time0, curve):
        if time0 is None:
            return curve
        m = curve.copy()
        m[t < time0] = np.nan
        return m

    # Find the earliest dose time across Vyvanse and Dex
    all_times = [td for td, _ in (VYVANSE or [])] + [td for td, _ in (DEX or [])]
    first_time = min(all_times) if all_times else None

    # PK curves (solid)
    total_pk_plot = mask_before(first_time, total_pk)
    total_pk_line, = ax.plot(t, total_pk_plot, linewidth=2.6, color=COLORS['total_pk'], label="Total (PK)")
    total_pk_color = total_pk_line.get_color()
    # Only plot Vyvanse PK if scheduled, and hide pre-dose baseline
    vyv_pk_line = None
    if VYVANSE:
        vyv_first = min(td for td, _ in VYVANSE)
        vyv_sum_plot = mask_before(vyv_first, vyv_sum)
        vyv_pk_line, = ax.plot(t, vyv_sum_plot, linewidth=2.0, color=COLORS['vyv_pk'], label="Vyvanse (PK)")
        vyv_pk_color = vyv_pk_line.get_color()
    else:
        vyv_pk_color = COLORS['vyv_pk']
    # Build effective Dex colors that avoid Total/Vyvanse colors
    reserved = {total_pk_color, vyv_pk_color}
    dex_colors = [c for c in DEX_BASE_COLORS if c not in reserved] or DEX_BASE_COLORS[:]
    dex_pk_colors, dex_pk_lines = [], []
    for i, ((td, d), curve) in enumerate(zip(DEX, dex_pk_curves)):
        col = dex_colors[i % len(dex_colors)]
        pk_line, = ax.plot(t, curve, linestyle="--", linewidth=1.5, color=col, label=f"Dex {d}mg @ {label_hour(td)} (PK)")
        dex_pk_lines.append(pk_line)
        dex_pk_colors.append(pk_line.get_color())

    # Perceived (PD) curves (dotted) — component-wise using kernel
    dt_h = t[1] - t[0]
    # Zero-DC biexponential kernel to produce faster wear-off perceived effect
    vyv_pd_components = [apply_pd_kernel(
                             np.nan_to_num(vc), dt_h,
                             PD.get('vyv_tau_r', 1.0), PD.get('vyv_tau_d', 6.0), PD.get('vyv_gain', 1.0),
                             match='peak', peak_scale=PD.get('pd_peak_scale', 1.0), clamp_scale=PD.get('pd_max_scale', None))
                         for (td, _), vc in zip(VYVANSE, vyv_pk_curves)] if vyv_pk_curves else []
    dex_pd_components = [apply_pd_kernel(
                             np.nan_to_num(dc), dt_h,
                             PD.get('dex_tau_r', 0.5), PD.get('dex_tau_d', 3.0), PD.get('dex_gain', 1.0),
                             match='peak', peak_scale=PD.get('pd_peak_scale', 1.0), clamp_scale=PD.get('pd_max_scale', None))
                         for (td, _), dc in zip(DEX, dex_pk_curves)]
    vyv_pd = np.nansum(np.vstack(vyv_pd_components), axis=0) if vyv_pd_components else np.zeros_like(t)
    total_pd = vyv_pd + (np.nansum(np.vstack(dex_pd_components), axis=0) if dex_pd_components else 0.0)

    # Apply visibility mask (keep Total perceived curve intact apart from the floor)
    floor = float(PD.get('pd_floor', 0.0))
    total_pd_m = mask_below(total_pd, floor)
    vyv_pd_m = mask_below(vyv_pd, floor)
    dex_pd_components_m = [mask_below(curve, floor) for curve in dex_pd_components]

    # Perceived component curves are plotted fully; no special masking is applied beyond the floor.
    # Mask perceived totals before first dose to avoid a zero baseline
    total_pd_plot = mask_before(first_time, total_pd_m)
    total_pd_line, = ax.plot(t, total_pd_plot, linewidth=2.6, linestyle=":", color=total_pk_color, alpha=0.9, label="Total (perceived)")
    vyv_pd_line = None
    if VYVANSE:
        vyv_pd_line, = ax.plot(t, vyv_pd_m, linewidth=2.0, linestyle=":", color=vyv_pk_color, alpha=0.85, label="Vyvanse (perceived)")
    dex_pd_lines = []
    for i, ((td, d), curve) in enumerate(zip(DEX, dex_pd_components_m)):
        col = dex_pk_colors[i] if i < len(dex_pk_colors) else (dex_colors[i % len(dex_colors)])
        pd_line, = ax.plot(t, curve, linestyle=":", linewidth=1.5, alpha=0.8, color=col, label=f"Dex {d}mg @ {label_hour(td)} (perceived)")
        dex_pd_lines.append(pd_line)

    # Stop-after projections
    def masked_from(time0, curve):
        m = curve.copy()
        m[t < time0] = np.nan
        return m

    # Do not render the last stop-after branch; it matches the final total
    stop_after_lines_pairs = []
    for i in range(max(0, len(DEX) - 1)):
        stop_pk = vyv_sum + sum(np.nan_to_num(c) for c in dex_pk_curves[:i+1])
        # PD branch: sum PD components up to i (plus Vyvanse PD)
        stop_pd = vyv_pd + np.nansum(np.vstack(dex_pd_components[:i+1]), axis=0)
        branch_time = DEX[i+1][0]
        col = dex_pk_colors[i] if i < len(dex_pk_colors) else (dex_colors[i % len(dex_colors)])
        sapk_line, = ax.plot(t, masked_from(branch_time, stop_pk), linestyle="--", linewidth=1.0, color=col, alpha=0.6, label=f"Stop after Dex {i+1} (PK)")
        sapd_line, = ax.plot(t, masked_from(branch_time, stop_pd), linestyle=":", linewidth=1.0, color=col, alpha=0.9, label=f"Stop after Dex {i+1} (perceived)")
        stop_after_lines_pairs.append((sapk_line, sapd_line))

    # Caffeine PK (separate right-hand y-axis; does not add to amphetamine totals)
    caffeine_total = None
    caffeine_line = None
    ax_caf = None
    if 'CAFFEINE' in globals() and CAFFEINE:
        caffeine_total = caffeine_total_curve(t, CAFFEINE)
        if caffeine_total is not None:
            ax_caf = ax.twinx()
            # De-emphasized but readable styling for caffeine scale/line
            caffeine_color = COLORS['caffeine']
            caffeine_line, = ax_caf.plot(t, caffeine_total, linewidth=1.9, linestyle='-', color=caffeine_color, alpha=0.75, label='Caffeine (PK)')
            # Configure caffeine axis on the LEFT (offset outward) alongside amphetamine scale
            caf_ymax = float(np.nanmax(caffeine_total)) if np.isfinite(np.nanmax(caffeine_total)) else 1.0
            ax_caf.set_ylim(0, np.ceil(caf_ymax * 1.1))
            # Hide right spine/ticks; show left spine offset outward to avoid overlap
            ax_caf.spines['right'].set_visible(False)
            ax_caf.spines['left'].set_visible(True)
            ax_caf.spines['left'].set_position(('outward', 42))
            ax_caf.spines['left'].set_color('#aaaaaa')
            ax_caf.yaxis.set_label_position('left')
            ax_caf.yaxis.tick_left()
            ax_caf.set_ylabel('Caffeine (mg, model)', color=caffeine_color)
            ax_caf.tick_params(axis='y', colors='#888888', labelsize=9)

    # Mark dose times (verticals matching dose colors) — skip first dose(s) of the day
    _EPS = 1e-6
    all_dose_times = [td for td, _ in (VYVANSE or [])] + [td for td, _ in (DEX or [])]
    first_time_global = min(all_dose_times) if all_dose_times else None
    for td, _ in VYVANSE:
        if first_time_global is not None and abs(td - first_time_global) < _EPS:
            continue
        ax.axvline(td, linestyle="--", linewidth=1.0, alpha=0.52, color=vyv_pk_color)
    for i, (td, _) in enumerate(DEX):
        if first_time_global is not None and abs(td - first_time_global) < _EPS:
            continue
        col = dex_pk_colors[i] if i < len(dex_pk_colors) else COLORS['neutral_marker']
        ax.axvline(td, linestyle="--", linewidth=1.0, alpha=0.52, color=col)

    # Axes/labels
    ax.set_xticks(range(int(t_start), int(t_end) + 1))
    ax.set_xticklabels([label_hour(h) for h in range(int(t_start), int(t_end) + 1)])
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_title(
        f"Vyvanse + Dex — PK (solid) vs perceived (dotted) | τr={PD.get('dex_tau_r',0.5)}h, τd={PD.get('dex_tau_d',3.0)}h, peak≈{PD.get('pd_peak_scale',1.0)}×PK, clamp≤{PD.get('pd_max_scale',1.1)}×PK"
    )
    ax.set_xlabel("Hour of Day")
    ax.set_ylabel("Amphetamine (PK/Perceived; a.u.)")
    ymax = float(np.nanmax(total_pk))
    ax.set_ylim(0, np.ceil(ymax * 1.1))
    ax.set_xlim(t_start, t_end)
    # Legend ordering: caffeine first, then amphetamine totals; then Vyvanse pair, Dex PK then Dex perceived; stop-after pairs at end
    handles = []
    if caffeine_line is not None:
        handles.append(caffeine_line)
    handles += [total_pk_line, total_pd_line]
    if vyv_pk_line is not None:
        handles.append(vyv_pk_line)
        if vyv_pd_line is not None:
            handles.append(vyv_pd_line)
    handles += dex_pk_lines + dex_pd_lines
    for sapk, sapd in stop_after_lines_pairs:
        handles += [sapk, sapd]
    ax.legend(handles=handles, labels=[h.get_label() for h in handles], fontsize=8, ncol=2)
    fig.tight_layout()
    return fig

# === Entrypoint ===
def run(save_fig=None):
    t_start, t_end = compute_time_window(VYVANSE, DEX, CAFFEINE)
    t = np.linspace(t_start, t_end, int((t_end - t_start) * 60 / RES_MIN) + 1)
    vyv_sum, vyv_pk_curves, dex_pk_curves, total_pk = build_pk(t)
    fig = plot_overlay(t, vyv_sum, vyv_pk_curves, dex_pk_curves, total_pk, default_PD, t_start, t_end)
    if save_fig:
        path = save_fig if save_fig != "default" else DEFAULT_OUTPUT
        save_figure_safely(fig, path)
    try:
        plt.show()
    except KeyboardInterrupt:
        plt.close('all')
        raise SystemExit(0)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        '--save-fig', nargs='?', const='default',
        help=f'Optionally save the chart to a file (format inferred from extension); omit path to use default: {DEFAULT_OUTPUT}'
    )
    args = p.parse_args()
    run(save_fig=args.save_fig)
