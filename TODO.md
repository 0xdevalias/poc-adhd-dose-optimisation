# TODO

- Refactor plotting scripts to DRY shared logic
  - Targets: `graph-dex-only-curves.py`, `graph-vyvanse-with-dex-curves.py`
  - Extract common helpers into a small module (e.g. `pk_plot_utils.py`):
    - `bateman`, `label_hour`, `curves_from_schedule`
    - Shared plotting style (grid, legends, dose marker styling)
  - Centralize pharmacokinetic config (Vyvanse + Dex) and `dex_mode_label` lines
    - Keep the simple toggle-by-comment pattern intact
  - Parameterize schedules/doses; expose functions that build curves and totals
  - Add CLI args

