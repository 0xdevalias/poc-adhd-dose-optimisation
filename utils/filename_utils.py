from __future__ import annotations

from datetime import date
from typing import Iterable, List, Sequence, Tuple, Union, Optional
from .dosing_utils import vyvanse_dex_eq_to_capsule_mg


Number = Union[int, float]
Schedule = Sequence[Tuple[Number, Number]]


def _fmt_num(x: Number) -> str:
    """Format a number compactly (trim trailing zeros, avoid sci notation)."""
    try:
        xf = float(x)
    except Exception:
        return str(x)
    # Use general format with sufficient precision, then trim trailing zeros/dot
    s = f"{xf:.10g}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


def build_schedule_filename(
    component: str,
    *,
    date_str: Optional[str] = None,
    vyvanse: Optional[Schedule] = None,
    dex: Optional[Schedule] = None,
    ext: str = "svg",
) -> str:
    """Construct a schedule-based filename like:

    YYYY-MM-DD-<component>-[vyv-<CAPmg>@<t>+...][<dexmg>@<t>+...].<ext>

    - `component` is a required descriptor (e.g., 'pk-vs-perceived').
    - Vyvanse entries labeled by inferred capsule mg (inverse of 0.4 ratio).
    - Dex entries labeled by their dose mg as provided.
    - Entries are ordered by time within each group; Vyvanse tokens come first if present.
    - If no doses are provided, returns: YYYY-MM-DD-<component>.<ext>
    """
    ds = date_str or date.today().isoformat()

    tokens: List[str] = []

    # Vyvanse tokens (use inferred capsule mg from dex-equivalent values)
    if vyvanse:
        vyv_sorted = sorted(vyvanse, key=lambda p: float(p[0]))
        for t_h, dex_eq_mg in vyv_sorted:
            cap_mg = vyvanse_dex_eq_to_capsule_mg(dex_eq_mg)
            tokens.append(f"vyv-{_fmt_num(cap_mg)}mg@{_fmt_num(t_h)}")

    # Dex tokens
    if dex:
        dex_sorted = sorted(dex, key=lambda p: float(p[0]))
        for t_h, mg in dex_sorted:
            tokens.append(f"{_fmt_num(mg)}mg@{_fmt_num(t_h)}")

    suffix = ("-" + "+".join(tokens)) if tokens else ""
    return f"{ds}-{component}{suffix}.{ext}"
