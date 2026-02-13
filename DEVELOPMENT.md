# CytoAtlas Report — Development Guide

## Structure

```
index.html                   Navigation landing page (GitHub Pages entry)
baseline/                    Original PI report (preserved)
  index.html                 16 MB interactive report with Plotly charts
  methodology.html           LinCytoSig methodology supplement
  REPORT.md / METHODS.md     Full report and methods documentation
reports/
  weekly/                    Quick weekly updates (~15 min to write)
  monthly/                   Synthesized from weekly reports
  quarterly/                 Major themes + manuscript readiness
  half-year/                 Near-manuscript-ready summaries
  manuscript/                Living manuscript draft + journal configs
figures/                     Shared figure library (15 baseline PNGs)
  weekly/                    New figures from weekly reports
  manuscript/                Curated manuscript figures
templates/                   Report templates (weekly, monthly, quarterly, half-year)
scripts/                     Build and generation scripts
_includes/                   Shared HTML/CSS for site generation
```

## Workflow

### Write a weekly report
```bash
python scripts/init_weekly.py              # current week
python scripts/init_weekly.py --week 2026-W08
```

### Synthesize higher-level reports
```bash
python scripts/synthesize.py monthly 2026-03
python scripts/synthesize.py quarterly 2026-Q2
python scripts/synthesize.py half-year 2026-H1
```

### Build HTML site
```bash
pip install markdown pyyaml jinja2
python scripts/build_site.py --all
python scripts/build_site.py --report reports/weekly/2026-W08.md
```

## Report Hierarchy

Weekly → Monthly → Quarterly → Half-Year → Manuscript

Each level synthesizes from the level below, with full traceability via source references. Templates include manuscript-mapping comments so content flows naturally into the final paper.

## Target Journals

Journal configurations in `reports/manuscript/journal-configs/`:
- Nature (3,500 words, 6 main figures)
- Nature Methods (3,000-5,000 words, Resource article type)
- Nature Medicine (4,000 words, Analysis article type)
- Cell (7,000 words, 7 figures, STAR Methods)
- Science (3,000 words, 3-5 figures)

## Baseline Documentation

| Document | Description |
|----------|-------------|
| [`baseline/REPORT.md`](baseline/REPORT.md) | Full baseline report |
| [`baseline/METHODS.md`](baseline/METHODS.md) | Per-figure provenance and statistical methods |
| [`baseline/DATA_DICTIONARY.md`](baseline/DATA_DICTIONARY.md) | Column definitions for all data files |
| [`baseline/PROVENANCE.md`](baseline/PROVENANCE.md) | Pipeline diagram and data flow |

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/generate_report_figures.py` | Generate 15 static matplotlib figures (existing) |
| `scripts/generate_interactive_report.py` | Generate interactive HTML report (existing) |
| `scripts/init_weekly.py` | Scaffold a new weekly report from template |
| `scripts/build_site.py` | Convert Markdown reports to HTML site |
| `scripts/synthesize.py` | Aggregate reports from lower to higher levels |
