#!/usr/bin/env python3

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

# === Dex-only scenario ===
t_dex      = [8.0,  9.5, 11.0, 11.75, 13.0]   # Dex dose times (h of day)
dex_mg     = [15.0, 5.0,  7.5,  2.5,   7.5]   # Dex dose amounts (mg)

dex_curves       = curves_from_schedule(t_dex, dex_mg, ka_ir, ke_ir)
total_dex_only   = sum(np.nan_to_num(c) for c in dex_curves)

# === Dynamic y-limit for full visibility with a little headroom ===
ymax = max(np.max(total_dex_only), np.max(total_ref))
y_top = np.ceil(ymax * 1.08)  # 8% headroom, rounded up

# === Plot ===
plt.figure(figsize=(13, 7))

# Reference total (dotted)
plt.plot(t, total_ref, linewidth=2.4, linestyle=":", label="Total (Vyvanse + Dex reference)")

# Dex-only components (dashed)
labels = [f"Dex-only {dose:g}mg @ {label_hour(td)}" for td, dose in zip(t_dex, dex_mg)]
dex_lines = []
for curve, lab in zip(dex_curves, labels):
    line, = plt.plot(t, curve, linestyle="--", linewidth=1.6, label=lab)
    dex_lines.append(line)

# Dex-only total (solid)
plt.plot(t, total_dex_only, linewidth=2.8, linestyle="-", label="Total (Dex-only model)")

# Dose markers for Dex-only (match line colors)
for td, line in zip(t_dex, dex_lines):
    plt.axvline(td, linestyle="--", linewidth=1.0, alpha=0.52, color=line.get_color())

# Hourly grid & ticks (8am → 8am next day)
xticks = list(range(8, 33, 1))
plt.xticks(xticks, [label_hour(h) for h in xticks], rotation=0)
plt.grid(True, which="both", axis="both", alpha=0.33, linestyle="--", linewidth=0.7)

plt.title(f"Dex-only Model vs Vyvanse+Dex Reference\nDex IR: {dex_mode_label}")
plt.xlabel("Hour of Day")
plt.ylabel("Relative Effect (arbitrary units)")
plt.ylim(0, y_top)
plt.xlim(8, 32)
plt.legend(ncol=2, fontsize=9)
plt.tight_layout()
plt.show()

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
