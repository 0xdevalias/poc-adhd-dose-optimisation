"""Centralized plotting style constants shared across scripts."""

# Base palette for per-dose Dex lines and related projections
DEX_BASE_COLORS = [
    "tab:purple", "tab:green", "tab:red",
    "tab:brown", "tab:pink", "tab:olive", "tab:cyan",
    "mediumpurple", "darkseagreen", "lightsalmon"
]

# Named colors used across plots
COLORS = {
    "total_pk": "tab:blue",     # also used by Total (perceived)
    "vyv_pk": "tab:orange",     # also used by Vyvanse (perceived)
    "neutral_marker": "tab:gray", # fallback for vertical markers
    "caffeine": "#555555",       # caffeine line color (de-emphasized dark gray)
}
