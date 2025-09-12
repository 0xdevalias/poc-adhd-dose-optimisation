#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import math
from utils.save_utils import save_figure_safely
from utils.pk_models import bateman, curves_from_schedule
from utils.plot_utils import format_time_12h
from utils.style import DEX_BASE_COLORS, COLORS
from utils.dosing_utils import vyvanse_cap_to_dex_eq
from utils.filename_utils import build_schedule_filename

def _build_schedule(times, doses):
    return list(zip(times, doses))

# === Pharmacokinetics ===
# (t½ = half-life)

# Vyvanse (LDX (lisdexamfetamine) → d-amphetamine) — reference
ka_vyv         = 0.80              # tuned for Tmax ≈ 3.5–4 h
ke_vyv         = np.log(2) / 11.0  # plasma half-life (t½) ≈ 10–12 h

# # Dex IR (actual PK; slower absorption / longer plasma half-life)
# ka_ir          = 1.0               # tuned for Tmax ≈ 2.5–3.5 h
# ke_ir          = np.log(2) / 11.0  # plasma half-life (t½) ≈ 10–12 h
# dex_mode_label = "PK (effective half-life) — ka=1.00, t1/2=11h"

# (Optional) Dex IR (perceived effect; faster absorption / shorter effective half-life)
ka_ir          = 1.4               # tuned for Tmax ≈ 1–2 h
ke_ir          = np.log(2) / 2.7   # effective half-life (t½) ≈ 3–4 h
dex_mode_label = "Perceived effect — ka=1.40, t1/2=2.7h"

# Output filename handling
DEFAULT_OUTPUT = "auto"  # sentinel; auto-generate from schedule unless overridden

# === Reference scenario: Vyvanse + Dex ===
# Use the same 'times' + 'doses' pattern for Vyvanse and Dex add-ons
t_vyv          = [8.0]                          # Vyvanse dose time(s) (h of day)
vyv_doses_mg   = [vyvanse_cap_to_dex_eq(30.0)]  # Vyvanse equivalent d-amphetamine (mg)

t_ref_dex      = [8.0, 11.0, 13.0]     # Dex top-up times (h of day)
ref_dex_mg     = [5.0, 5.0, 5.0]       # Dex top-up doses (mg)

# === Dex-only scenario ===
t_dex      = [8.0,  9.5, 11.0, 11.75, 13.0]   # Dex dose times (h of day)
dex_mg     = [15.0, 5.0,  7.5,  2.5,   7.5]   # Dex dose amounts (mg)

# === Dynamic 24h window based on earliest dose across all plotted schedules ===
all_times = []
all_times += list(t_vyv or [])
all_times += list(t_ref_dex or [])
all_times += list(t_dex or [])
if all_times:
    start_h = float(math.floor(min(all_times)))
else:
    start_h = 8.0
end_h = start_h + 24.0
t = np.linspace(start_h, end_h, int((end_h - start_h) * 60) + 1)  # 1-min resolution
vyv_curves_ref = curves_from_schedule(t, _build_schedule(t_vyv, vyv_doses_mg), ka_vyv, ke_vyv)
ref_dex_curves = curves_from_schedule(t, _build_schedule(t_ref_dex, ref_dex_mg), ka_ir,  ke_ir)

vyv_curve_ref  = sum(np.nan_to_num(v) for v in vyv_curves_ref)
total_ref      = vyv_curve_ref + sum(np.nan_to_num(c) for c in ref_dex_curves)

dex_curves       = curves_from_schedule(t, _build_schedule(t_dex, dex_mg), ka_ir, ke_ir)
total_dex_only   = sum(np.nan_to_num(c) for c in dex_curves)

# Mask helper to hide pre-dose zero baselines for plotted totals
def mask_before(time0, curve):
    m = curve.copy()
    m[t < time0] = np.nan
    return m

first_time = min(all_times) if all_times else None
total_ref_plot = mask_before(first_time, total_ref) if first_time is not None else total_ref
total_dex_only_plot = mask_before(first_time, total_dex_only) if first_time is not None else total_dex_only

# Helper to mask a curve so it only shows from a given time onwards
def mask_from(time0, curve):
    m = curve.copy()
    m[t < time0] = np.nan
    return m

# === Dynamic y-limit for full visibility with a little headroom ===
ymax = max(np.max(total_dex_only), np.max(total_ref))
y_top = np.ceil(ymax * 1.08)  # 8% headroom, rounded up

# === Plot ===
fig, ax = plt.subplots(figsize=(13, 7))

# Reference total (dotted)
ref_total_line, = ax.plot(t, total_ref_plot, linewidth=2.4, linestyle=":", color=COLORS["total_pk"], label="Total (Vyvanse + Dex reference)")

# Dex-only components (dashed)
labels = [f"Dex-only {dose:g}mg @ {format_time_12h(td)}" for td, dose in zip(t_dex, dex_mg)]
dex_lines = []
for i, (curve, lab) in enumerate(zip(dex_curves, labels)):
    col = DEX_BASE_COLORS[i % len(DEX_BASE_COLORS)]
    line, = ax.plot(t, curve, linestyle="--", linewidth=1.6, color=col, label=lab)
    dex_lines.append(line)

# Stop-after projections (dotted, behind totals) — one per Dex dose except last
# Skip the first branch in Dex-only, as it duplicates the first-dose curve.
stop_after_lines = []
for i in range(max(0, len(dex_curves) - 1)):
    if i == 0:
        continue
    branch_time = t_dex[i + 1]
    included_time = t_dex[i]
    partial = sum(np.nan_to_num(c) for c in dex_curves[: i + 1])
    color = dex_lines[i].get_color() if i < len(dex_lines) else None
    line, = ax.plot(
        t,
        mask_from(branch_time, partial),
        linestyle=":",
        linewidth=1.4,
        alpha=0.85,
        color=color,
        label=f"Stop after {format_time_12h(included_time)} Dex"
    )
    stop_after_lines.append(line)

# Dex-only total (solid)
dex_total_line, = ax.plot(t, total_dex_only_plot, linewidth=2.8, linestyle="-", color=COLORS["total_pk"], label="Total (Dex-only model)")

# Dose markers for Dex-only (match line colors) — skip first Dex dose of the day
_EPS = 1e-6
first_dex_time = min(t_dex) if t_dex else None
for td, line in zip(t_dex, dex_lines):
    if first_dex_time is not None and abs(td - first_dex_time) < _EPS:
        continue
    ax.axvline(td, linestyle="--", linewidth=1.0, alpha=0.52, color=line.get_color())

# Hourly grid & ticks (start → start next day)
xticks = list(range(int(start_h), int(end_h) + 1, 1))
ax.set_xticks(xticks)
ax.set_xticklabels([format_time_12h(h) for h in xticks], rotation=0)
ax.grid(True, which="both", axis="both", alpha=0.33, linestyle="--", linewidth=0.7)

ax.set_title(f"Dex-only Model vs Vyvanse+Dex Reference\nDex IR: {dex_mode_label}")
ax.set_xlabel("Hour of Day")
ax.set_ylabel("Relative Effect (arbitrary units)")
ax.set_ylim(0, y_top)
ax.set_xlim(start_h, end_h)
# Legend ordering: place stop-after entries at end of legend
handles = [ref_total_line, dex_total_line] + dex_lines + stop_after_lines
ax.legend(handles=handles, labels=[h.get_label() for h in handles], ncol=2, fontsize=9)
fig.tight_layout()

# Auto-generate DEFAULT_OUTPUT from schedule only if left as the 'auto' sentinel
if isinstance(DEFAULT_OUTPUT, str) and DEFAULT_OUTPUT == "auto":
    DEFAULT_OUTPUT = build_schedule_filename(
        'dex-only-curves',
        vyvanse=_build_schedule(t_vyv, vyv_doses_mg),
        dex=_build_schedule(t_dex, dex_mg),
        ext="svg",
    )

# Optional: save the chart to a file (format inferred from extension)
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--save-fig",
    nargs="?",
    const="default",
    metavar="PATH",
    help="Optionally save the chart to a file. Provide a path or omit for a default filename."
)
try:
    args, _ = parser.parse_known_args()
    save_path = None
    if getattr(args, "save_fig", None) is not None:
        save_path = args.save_fig if args.save_fig != "default" else DEFAULT_OUTPUT
    if save_path:
        save_figure_safely(fig, save_path)
except SystemExit:
    pass  # ignore argparse errors in interactive contexts

try:
    plt.show()
except KeyboardInterrupt:
    plt.close('all')
    raise SystemExit(0)

# === Values at key targets for Dex-only ===
def value_at_hour(hour, curve):
    """Interpolate value at given hour from a single total curve."""
    return float(np.interp(hour, t, curve))

vals = {
    "10am": value_at_hour(10.0, total_dex_only),
    "12pm": value_at_hour(12.0, total_dex_only),
    "2pm":  value_at_hour(14.0, total_dex_only),
}

print("Dex-only model totals at targets:")
for k, v in vals.items():
    print(f"  {k}: {v:.2f}")
