from typing import List, Tuple, Optional


def expand_split_dose(start_time_h: float, total_amount: float, duration_min: float, *, parts: Optional[int] = None) -> List[Tuple[float, float]]:
    """Split a dose across a duration (minutes), optionally into a fixed number of parts.

    - start_time_h: when ingestion begins (hours)
    - total_amount: total dose amount (same units you model, e.g., mg)
    - duration_min: how long ingestion takes (minutes). If <= 0, returns a single dose.
    - parts: optional explicit number of equal sub-doses across `duration_min`.
             If not provided, defaults to one sub-dose per minute.

    Returns list of (time_h, amount) pairs summing to total_amount.
    """
    if duration_min is None or duration_min <= 0:
        return [(float(start_time_h), float(total_amount))]
    # Determine number of sub-doses
    n_parts = int(parts) if parts is not None else int(round(float(duration_min)))
    n_parts = max(1, n_parts)
    step_h = (float(duration_min) / n_parts) / 60.0
    amount_each = float(total_amount) / n_parts
    return [(float(start_time_h) + i * step_h, amount_each) for i in range(n_parts)]


def shots_to_caffeine_mg(shots: float, *, mg_per_shot: float = 75.0) -> float:
    """Estimate caffeine (mg) from espresso shots.

    Default is 75 mg per single shot (midpoint of common 60–80 mg range).
    Adjust `mg_per_shot` to match your cafe/coffee setup.

    Reference ranges (per serving):
      - Espresso (single, 30 mL): 60–80 mg
      - Espresso (double, 60 mL): 120–160 mg
      - Lungo (60–90 mL): 80–100 mg
      - Ristretto (≈20 mL): 50–60 mg
      - Drip/Filter (240 mL): 95–200 mg
      - French press (240 mL): 80–135 mg
      - Instant (240 mL): 60–90 mg
      - Cold brew concentrate (240 mL): 150–300 mg (can be higher)
      - Decaf (240 mL): 2–15 mg
    """
    return float(shots) * float(mg_per_shot)


def aeropress_scoops_to_caffeine_mg(
    scoops: float,
    *,
    grams_per_scoop: float = 13.0,
    mg_per_gram: float = 12.0,
) -> float:
    """Estimate caffeine (mg) from AeroPress scoops via bean weight.

    - grams_per_scoop: default 13 g (AeroPress scoop ≈ 11–14 g; bean size/roast vary).
    - mg_per_gram: default 12 mg/g (Arabica ~10–12 mg/g; Robusta ~18–20 mg/g).

    Note: AeroPress extraction is highly variable; tune both parameters for your recipe/beans.

    Reference (AeroPress scoop and espresso comparisons):
      - Scoop volume: ≈ 11–14 g whole beans (depends on bean size, roast, grind).
      - Light, dense beans → closer to 14 g; dark, puffier beans → closer to 11 g.
      - Roughly ≈ 1.5 standard tablespoons by volume.
      - Espresso beans per shot: single ≈ 7–9 g; double ≈ 14–18 g.
      - One AeroPress scoop (≈ 11–14 g) ≈ a strong single or lighter double by bean mass.
      - Practical: 1–1.5 scoops for a strong mug; ~2 scoops for a large mug or to split.
    """
    grams = float(scoops) * float(grams_per_scoop)
    return grams_to_caffeine_mg(grams, mg_per_gram=mg_per_gram)


def grams_to_caffeine_mg(grams: float, *, mg_per_gram: float = 12.0) -> float:
    """Estimate caffeine (mg) from grams of dry coffee beans.

    - mg_per_gram sets the assumed caffeine yield per gram of beans.
      Typical rules of thumb: Arabica ≈ 10–12 mg/g, Robusta ≈ 18–20 mg/g.
    - Default 12 mg/g is a pragmatic Arabica midpoint; tune to your beans/recipe.
    - Extraction varies with brew method, grind, roast, water ratio, and time —
      treat this as a coarse estimate and calibrate `mg_per_gram` if you have measurements.
    - Quick examples (Arabica): 7 g → ~70–85 mg; 14 g → ~140–170 mg; 13 g (≈ 1 AeroPress scoop) → ~130–155 mg.

    Reference (bean measures):
      - Espresso beans per shot: single ≈ 7–9 g; double ≈ 14–18 g.
      - AeroPress scoop: ≈ 11–14 g whole beans (bean size/roast/grind vary).
    """
    return float(grams) * float(mg_per_gram)
