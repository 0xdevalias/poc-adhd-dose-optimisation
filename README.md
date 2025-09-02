# poc-adhd-dose-optimisation

Proof‑of‑concept (PoC) scripts for exploring ADHD medication dosing patterns and visualising relative effect over time. These models are intentionally simple and tuned for quick iteration, not scientific accuracy.

> [!IMPORTANT]
> The dose impact and decay models here may not match real pharmacokinetics or perceived effects. Use as exploratory visualisations only — not medical advice.

## Notes and Limitations

- Simplified one‑compartment Bateman model; parameters are hand‑tuned for plausible shapes.
- “Perceived effect” uses a shorter effective half‑life to better match how IR Dex often feels, not measured plasma kinetics.
- Doses and schedules in scripts are examples only.
- Visual styles include a subtle dashed grid and faint color‑matched vertical lines at Dex dose times to aid reading without distraction.

### Toggling Dex IR Models (Perceived vs PK)

Both scripts include two Dex IR model blocks in the “Pharmacokinetics” section. Toggle by commenting/uncommenting the block you want. The active mode appears in the chart title.

Actual PK (longer plasma half‑life):

```python
# # Dex IR (actual PK; slower absorption / longer plasma half-life)
# ka_ir          = 1.0               # tuned for Tmax ≈ 2.5–3.5 h
# ke_ir          = np.log(2) / 11.0  # plasma half-life (t½) ≈ 10–12 h
# dex_mode_label = "PK (effective half-life) — ka=1.00, t1/2=11h"
```

Perceived effect (faster onset, shorter felt duration):

```python
# (Optional) Dex IR (perceived effect; faster absorption / shorter effective half-life)
ka_ir          = 1.4               # tuned for Tmax ≈ 1–2 h
ke_ir          = np.log(2) / 2.7   # effective half-life (t½) ≈ 3–4 h
dex_mode_label = "Perceived effect — ka=1.40, t1/2=2.7h"
```

## What’s Included

- `graph-vyvanse-with-dex-curves.py`: Models Vyvanse with Dex IR top‑ups (the “reference” scenario) with stop‑after projections.
- `graph-dex-only-curves.py`: Models a Dex IR‑only schedule and compares it against a reference Vyvanse+Dex total.

Both plotting scripts share the same structure and section order to keep diffs clean: helpers → time axis → pharmacokinetics → reference/schedule → y‑limits → plot → values at key targets.

## Environment Setup

You can use either `pyenv` + `pyenv-virtualenv` (preferred) or a standard Python `venv`.

### Option A: `pyenv` + `pyenv-virtualenv` (preferred)

If needed, install via Homebrew on macOS (https://brew.sh/):

```sh
brew update

brew install pyenv pyenv-virtualenv

# Add to your shell profile (zsh shown):
echo 'eval "$(pyenv init -)"' >> ~/.zprofile
echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.zprofile
exec $SHELL -l
```

1) Choose an installed Python version (or install a new one):

```sh
pyenv versions        # list installed
pyenv install -l      # list all available (if you need another)
```

2) Create a project virtualenv and set it locally:

```sh
pyenv virtualenv <python-version> poc-adhd-dose-optimisation
pyenv local poc-adhd-dose-optimisation   # writes .python-version
```

3) Install dependencies:

```sh
python -m pip install -r requirements.txt
```

### Option B: Standard virtualenv

```sh
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

## Running the Scripts

```sh
python graph-dex-only-curves.py
python graph-vyvanse-with-dex-curves.py
```

Each script opens a Matplotlib window with the plotted curves and prints values at key target times in the terminal.

## License

See the [LICENSE](LICENSE) file in this repository.
