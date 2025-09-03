#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

# === Helpers ===
def bateman(t, dose, t0, ka, ke):
    """One-compartment Bateman function for an oral dose starting at t0."""
    c = np.zeros_like(t)
    m = t >= t0
    tm = t[m] - t0
    if abs(ka - ke) < 1e-6:  # avoid division by ~0 in ka/(ka-ke)
        ka += 1e-6
    c[m] = dose * (ka / (ka - ke)) * (np.exp(-ke * tm) - np.exp(-ka * tm))
    c[c < 0] = 0
    return c

def curves_from_schedule(times, doses, ka, ke):
    """Return list of Bateman curves (one per (time,dose)) with pre-dose masked."""
    curves = [bateman(t, d, td, ka, ke) for td, d in zip(times, doses)]
    for curve, td in zip(curves, times):
        curve[t < td] = np.nan  # hide pre-dose segment for clarity
    return curves

def label_hour(h):
    """Format a decimal hour (0–24+) as a 12-hour label like '8am' or '1pm'."""
    h24 = int(h % 24)
    suffix = "am" if h24 < 12 else "pm"
    h12 = h24 % 12
    if h12 == 0:
        h12 = 12
    return f"{h12}{suffix}"

# === Time axis: full 24h from 08:00 today to 08:00 next day ===
t = np.linspace(8, 32, 24*60+1)  # 1-min resolution

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

# === Reference scenario: Vyvanse + Dex ===
# Use the same 'times' + 'doses' pattern for Vyvanse and Dex add-ons
t_vyv          = [8.0]                 # Vyvanse dose time(s) (h of day)
vyv_doses_mg   = [12.0]                # Vyvanse equivalent d-amphetamine (mg)

t_ref_dex      = [8.0, 11.0, 13.0]     # Dex top-up times (h of day)
ref_dex_mg     = [5.0, 5.0, 5.0]       # Dex top-up doses (mg)

vyv_curves_ref = curves_from_schedule(t_vyv,     vyv_doses_mg, ka_vyv, ke_vyv)
ref_dex_curves = curves_from_schedule(t_ref_dex, ref_dex_mg,   ka_ir,  ke_ir)

vyv_curve_ref  = sum(np.nan_to_num(v) for v in vyv_curves_ref)
total_ref      = vyv_curve_ref + sum(np.nan_to_num(c) for c in ref_dex_curves)

# Stop-after projections (drawn behind, color-matched to dose)
total_stop_after_dex1 = vyv_curve_ref + np.nan_to_num(ref_dex_curves[0])
total_stop_after_dex2 = vyv_curve_ref + np.nan_to_num(ref_dex_curves[0]) + np.nan_to_num(ref_dex_curves[1])

# === Dynamic y-limit for full visibility with a little headroom ===
ymax = np.max(total_ref)
y_top = np.ceil(ymax * 1.08)  # 8% headroom, rounded up

# === Plot ===
plt.figure(figsize=(13, 7))

# Projections first (so total overlays overlaps)
plt.plot(t, total_stop_after_dex1, linestyle=":", linewidth=1.8, color="tab:purple", alpha=0.85, label="Stop after 08:00 Dex")
plt.plot(t, total_stop_after_dex2, linestyle=":", linewidth=1.8, color="tab:green",  alpha=0.85, label="Stop after 11:00 Dex")

# Reference total (solid)
plt.plot(t, total_ref,  linewidth=2.6, color="tab:blue", label="Total (Vyvanse + Dex)")

# Base components (dashed)
plt.plot(t, vyv_curve_ref,  linewidth=2.0, color="orange", label="Vyvanse 30mg → dex (eq. 12mg) [Tmax≈3.5–4h]")
labels = [f"Dex IR {dose:g}mg @ {label_hour(td)}" for td, dose in zip(t_ref_dex, ref_dex_mg)]
colors = ["tab:purple", "tab:green", "gold"]
dex_lines = []
for curve, lab, col in zip(ref_dex_curves, labels, colors):
    line, = plt.plot(t, curve, linestyle="--", linewidth=1.8, color=col, label=f"{lab} ({'perceived' if 'Perceived' in dex_mode_label else 'PK'})")
    dex_lines.append(line)

# Dose markers (Dex only; match line colors)
for td, line in zip(t_ref_dex, dex_lines):
    plt.axvline(td, linestyle="--", linewidth=1.0, alpha=0.52, color=line.get_color())

# Hourly grid & ticks (8am → 8am next day)
xticks = list(range(8, 33, 1))
plt.xticks(xticks, [label_hour(h) for h in xticks], rotation=0)
plt.grid(True, which="both", axis="both", alpha=0.33, linestyle="--", linewidth=0.7)

plt.title(f"Vyvanse + Dex Model\nDex IR: {dex_mode_label}")
plt.xlabel("Hour of Day")
plt.ylabel("Relative Effect (arbitrary units)")
plt.ylim(0, y_top)
plt.xlim(8, 32)
plt.legend(ncol=2, fontsize=9)
plt.tight_layout()

# Optional: save an SVG of the chart
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--save-svg", nargs="?", const="default", metavar="PATH",
                    help="Optionally save the chart as an SVG. Provide a path or omit for a default filename.")
try:
    args, _ = parser.parse_known_args()
    svg_path = None
    if getattr(args, "save_svg", None) is not None:
        svg_path = args.save_svg if args.save_svg != "default" else "graph-vyvanse-with-dex-curves.svg"
    if svg_path:
        plt.savefig(svg_path, format="svg", bbox_inches="tight")
except SystemExit:
    pass  # ignore argparse errors in interactive contexts

plt.show()

# === Values at key targets for Vyvanse+Dex ===
def value_at_hour(hour, curve):
    """Interpolate value at given hour from a single total curve."""
    return float(np.interp(hour, t, curve))

vals = {
    "10am": value_at_hour(10.0, total_ref),
    "12pm": value_at_hour(12.0, total_ref),
    "2pm":  value_at_hour(14.0, total_ref),
}

print("Vyvanse+Dex totals at targets:")
for k, v in vals.items():
    print(f"  {k}: {v:.2f}")
