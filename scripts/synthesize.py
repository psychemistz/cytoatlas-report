#!/usr/bin/env python3
"""Semi-automated report synthesis: aggregate lower-level reports into
higher-level templates.

Does mechanical aggregation (section extraction, grouping by theme);
the researcher does intellectual synthesis by editing the output.

Usage:
    python scripts/synthesize.py monthly 2026-03
    python scripts/synthesize.py quarterly 2026-Q2
    python scripts/synthesize.py half-year 2026-H1
    python scripts/synthesize.py manuscript

Dependencies: pip install pyyaml
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    print("Missing dependency: pyyaml")
    print("Install with: pip install pyyaml")
    raise SystemExit(1)

ROOT = Path(__file__).resolve().parent.parent
REPORTS_DIR = ROOT / "reports"
TEMPLATES_DIR = ROOT / "templates"

# Project vocabulary for theme grouping
THEME_KEYWORDS = {
    "validation": [
        "validation", "correlat", "spearman", "rho", "benchmark",
        "ground truth",
    ],
    "cross-atlas": [
        "cross-atlas", "consistency", "atlas", "generalizab",
    ],
    "method-comparison": [
        "cytosig", "lincytosig", "secact", "method comparison",
        "10-way", "advantage",
    ],
    "cell-type": [
        "cell type", "cell-type", "celltype", "specificity",
        "t cell", "monocyte", "macrophage", "nk cell",
    ],
    "bulk-validation": [
        "bulk", "gtex", "tcga", "rna-seq",
    ],
    "biological": [
        "biological", "pathway", "cytokine", "interferon", "tnf",
        "il-", "tgfb", "ifn",
    ],
}


def parse_frontmatter(text: str) -> tuple[dict, str]:
    if text.startswith("---"):
        parts = text.split("---", 2)
        if len(parts) >= 3:
            meta = yaml.safe_load(parts[1]) or {}
            body = parts[2]
            return meta, body
    return {}, text


def extract_sections(text: str) -> dict[str, str]:
    """Extract content under each ## heading."""
    sections: dict[str, str] = {}
    current = None
    lines: list[str] = []

    for line in text.splitlines():
        if line.startswith("## "):
            if current is not None:
                sections[current] = "\n".join(lines).strip()
            current = line[3:].strip()
            lines = []
        else:
            lines.append(line)

    if current is not None:
        sections[current] = "\n".join(lines).strip()

    return sections


def classify_theme(text: str) -> list[str]:
    """Match text against project vocabulary, return matching themes."""
    text_lower = text.lower()
    matches = []
    for theme, keywords in THEME_KEYWORDS.items():
        if any(kw in text_lower for kw in keywords):
            matches.append(theme)
    return matches or ["general"]


def get_weeks_in_month(month_str: str) -> list[Path]:
    """Find weekly reports that fall within a given month (YYYY-MM)."""
    year, month = int(month_str[:4]), int(month_str[5:7])
    weekly_dir = REPORTS_DIR / "weekly"
    if not weekly_dir.exists():
        return []

    results = []
    for f in sorted(weekly_dir.glob("*.md")):
        raw = f.read_text(encoding="utf-8")
        meta, _ = parse_frontmatter(raw)
        date_range = meta.get("date_range", "")
        # Check if the week's Monday falls in this month
        if date_range:
            start_date = date_range.split(" to ")[0].strip()
            try:
                dt = datetime.strptime(start_date, "%Y-%m-%d")
                if dt.year == year and dt.month == month:
                    results.append(f)
            except ValueError:
                pass
        else:
            # Fallback: parse from filename YYYY-WXX
            week_str = f.stem
            if week_str.startswith(str(year)):
                results.append(f)

    return results


def get_months_in_quarter(quarter_str: str) -> list[Path]:
    """Find monthly reports for a given quarter (YYYY-QX)."""
    year = int(quarter_str[:4])
    q = int(quarter_str[-1])
    months = [(q - 1) * 3 + i for i in range(1, 4)]

    monthly_dir = REPORTS_DIR / "monthly"
    if not monthly_dir.exists():
        return []

    results = []
    for m in months:
        p = monthly_dir / f"{year}-{m:02d}.md"
        if p.exists():
            results.append(p)
    return results


def get_quarters_in_half(half_str: str) -> list[Path]:
    """Find quarterly reports for a given half-year (YYYY-HX)."""
    year = int(half_str[:4])
    h = int(half_str[-1])
    qs = [1, 2] if h == 1 else [3, 4]

    quarterly_dir = REPORTS_DIR / "quarterly"
    if not quarterly_dir.exists():
        return []

    results = []
    for q in qs:
        p = quarterly_dir / f"{year}-Q{q}.md"
        if p.exists():
            results.append(p)
    return results


def synthesize_monthly(month_str: str):
    """Generate a pre-filled monthly template from weekly reports."""
    weeks = get_weeks_in_month(month_str)
    if not weeks:
        print(f"No weekly reports found for {month_str}")
        print("Creating empty monthly template anyway...")

    template = (TEMPLATES_DIR / "monthly.md").read_text(encoding="utf-8")
    year, month = int(month_str[:4]), int(month_str[5:7])

    # Collect all weekly content
    all_findings: list[dict] = []
    all_methods: list[str] = []
    all_figures: list[dict] = []
    week_ids: list[str] = []

    for wf in weeks:
        raw = wf.read_text(encoding="utf-8")
        meta, body = parse_frontmatter(raw)
        sections = extract_sections(body)
        week_id = meta.get("week", wf.stem)
        week_ids.append(week_id)

        # Key Results
        if "Key Results" in sections:
            themes = classify_theme(sections["Key Results"])
            all_findings.append({
                "week": week_id,
                "content": sections["Key Results"],
                "themes": themes,
            })

        # Methods
        if "Methods & Technical" in sections:
            content = sections["Methods & Technical"].strip()
            if content and content != "-":
                all_methods.append(f"<!-- Source: {week_id} -->\n{content}")

        # Figures
        if "Figures" in sections:
            all_figures.append({
                "week": week_id,
                "content": sections["Figures"],
            })

    # Fill template
    output = template
    output = output.replace("YYYY-MM", month_str)
    output = output.replace(
        "YYYY-MM-DD to YYYY-MM-DD",
        f"{year}-{month:02d}-01 to {year}-{month:02d}-28",
    )
    output = output.replace("{{ month }}", month_str)
    output = output.replace(
        "{{ date_range }}",
        f"{year}-{month:02d}-01 to {year}-{month:02d}-28",
    )
    output = output.replace(
        '{{ source_weeks | join(", ") }}',
        ", ".join(week_ids),
    )
    output = output.replace("source_weeks: []", f"source_weeks: {week_ids}")

    # Group findings by theme
    if all_findings:
        themed: dict[str, list[dict]] = {}
        for finding in all_findings:
            for theme in finding["themes"]:
                themed.setdefault(theme, []).append(finding)

        findings_text = ""
        for i, (theme, items) in enumerate(themed.items(), 1):
            findings_text += f"### Theme {i}: {theme.replace('-', ' ').title()}\n"
            for item in items:
                findings_text += f"<!-- Source: {item['week']} -->\n"
                findings_text += item["content"] + "\n\n"

        output = output.replace(
            "### Theme 1: [title]\n<!-- Source: W?? -->\n\n-\n\n"
            "### Theme 2: [title]\n<!-- Source: W?? -->\n\n-",
            findings_text.strip(),
        )

    # Methods
    if all_methods:
        methods_text = "\n\n".join(all_methods)
        output = output.replace(
            "<!-- Deduplicated from weekly Methods & Technical sections -->\n\n-",
            f"<!-- Deduplicated from weekly Methods & Technical sections -->\n\n{methods_text}",
        )

    out_path = REPORTS_DIR / "monthly" / f"{month_str}.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(output, encoding="utf-8")
    print(f"Created: {out_path}")
    print(f"Source weeks: {week_ids}")
    print(">>> Edit the generated file to add intellectual synthesis <<<")


def synthesize_quarterly(quarter_str: str):
    """Generate a pre-filled quarterly template from monthly reports."""
    months = get_months_in_quarter(quarter_str)
    if not months:
        print(f"No monthly reports found for {quarter_str}")
        print("Creating empty quarterly template anyway...")

    template = (TEMPLATES_DIR / "quarterly.md").read_text(encoding="utf-8")
    year = int(quarter_str[:4])
    q = int(quarter_str[-1])
    month_ids = [m.stem for m in months]

    output = template
    output = output.replace("YYYY-QX", quarter_str)
    output = output.replace("{{ quarter }}", quarter_str)
    output = output.replace(
        '{{ source_months | join(", ") }}',
        ", ".join(month_ids),
    )
    output = output.replace("source_months: []", f"source_months: {month_ids}")

    # Extract key findings from monthlies
    if months:
        themes_text = ""
        theme_num = 1
        for mf in months:
            raw = mf.read_text(encoding="utf-8")
            meta, body = parse_frontmatter(raw)
            sections = extract_sections(body)
            if "Key Findings" in sections:
                themes_text += f"### {theme_num}. From {mf.stem}\n"
                themes_text += f"<!-- Source: {mf.stem} -->\n"
                themes_text += sections["Key Findings"] + "\n\n"
                theme_num += 1

        if themes_text:
            output = output.replace(
                "### 1. [Theme title]\n<!-- Synthesized from monthly key findings -->\n\n\n\n"
                "### 2. [Theme title]\n\n\n\n### 3. [Theme title]",
                themes_text.strip(),
            )

    start_month = (q - 1) * 3 + 1
    end_month = q * 3
    output = output.replace(
        "YYYY-MM-DD to YYYY-MM-DD",
        f"{year}-{start_month:02d}-01 to {year}-{end_month:02d}-28",
    )
    output = output.replace(
        "{{ date_range }}",
        f"{year}-{start_month:02d}-01 to {year}-{end_month:02d}-28",
    )

    out_path = REPORTS_DIR / "quarterly" / f"{quarter_str}.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(output, encoding="utf-8")
    print(f"Created: {out_path}")
    print(f"Source months: {month_ids}")
    print(">>> Edit the generated file to add intellectual synthesis <<<")


def synthesize_half_year(half_str: str):
    """Generate a pre-filled half-year template from quarterly reports."""
    quarters = get_quarters_in_half(half_str)
    if not quarters:
        print(f"No quarterly reports found for {half_str}")
        print("Creating empty half-year template anyway...")

    template = (TEMPLATES_DIR / "half-year.md").read_text(encoding="utf-8")
    year = int(half_str[:4])
    h = int(half_str[-1])
    quarter_ids = [q.stem for q in quarters]

    output = template
    output = output.replace("YYYY-HX", half_str)
    output = output.replace("{{ half_year }}", half_str)
    output = output.replace(
        '{{ source_quarters | join(", ") }}',
        ", ".join(quarter_ids),
    )
    output = output.replace(
        "source_quarters: []", f"source_quarters: {quarter_ids}"
    )

    if h == 1:
        output = output.replace(
            "YYYY-MM-DD to YYYY-MM-DD",
            f"{year}-01-01 to {year}-06-30",
        )
        output = output.replace(
            "{{ date_range }}",
            f"{year}-01-01 to {year}-06-30",
        )
    else:
        output = output.replace(
            "YYYY-MM-DD to YYYY-MM-DD",
            f"{year}-07-01 to {year}-12-31",
        )
        output = output.replace(
            "{{ date_range }}",
            f"{year}-07-01 to {year}-12-31",
        )

    # Pull quarter summaries
    if quarters:
        results_text = ""
        for qf in quarters:
            raw = qf.read_text(encoding="utf-8")
            meta, body = parse_frontmatter(raw)
            sections = extract_sections(body)
            if "Quarter Summary" in sections:
                results_text += f"### From {qf.stem}\n"
                results_text += sections["Quarter Summary"] + "\n\n"

        if results_text:
            output = output.replace(
                "### [Result theme 1]\n<!-- Organized as manuscript paragraphs -->\n\n\n\n"
                "### [Result theme 2]\n\n\n\n### [Result theme 3]",
                results_text.strip(),
            )

    out_path = REPORTS_DIR / "half-year" / f"{half_str}.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(output, encoding="utf-8")
    print(f"Created: {out_path}")
    print(f"Source quarters: {quarter_ids}")
    print(">>> Edit the generated file to add intellectual synthesis <<<")


def synthesize_manuscript():
    """Generate manuscript sections by pulling from the most recent
    half-year or quarterly reports."""
    # Find most recent reports at each level
    half_year_dir = REPORTS_DIR / "half-year"
    quarterly_dir = REPORTS_DIR / "quarterly"
    weekly_dir = REPORTS_DIR / "weekly"

    sources = []
    for d in [half_year_dir, quarterly_dir, weekly_dir]:
        if d.exists():
            sources.extend(sorted(d.glob("*.md"), reverse=True))

    if not sources:
        print("No reports found to synthesize into manuscript.")
        print("Write some weekly reports first!")
        return

    draft_path = REPORTS_DIR / "manuscript" / "draft.md"
    if draft_path.exists():
        print(f"Manuscript draft already exists: {draft_path}")
        print("Not overwriting. Edit it directly or remove it first.")
        return

    print(f"Found {len(sources)} source reports.")
    print("Manuscript scaffold should be seeded from baseline content.")
    print(f"See: {draft_path}")
    print(">>> Use the baseline report sections to populate the draft <<<")


def main():
    parser = argparse.ArgumentParser(
        description="Synthesize reports from lower to higher levels"
    )
    parser.add_argument(
        "level",
        choices=["monthly", "quarterly", "half-year", "manuscript"],
        help="Target report level",
    )
    parser.add_argument(
        "period",
        nargs="?",
        help="Period identifier (e.g. 2026-03, 2026-Q2, 2026-H1)",
    )
    args = parser.parse_args()

    if args.level == "manuscript":
        synthesize_manuscript()
    elif args.period is None:
        parser.error(f"{args.level} requires a period argument")
    elif args.level == "monthly":
        synthesize_monthly(args.period)
    elif args.level == "quarterly":
        synthesize_quarterly(args.period)
    elif args.level == "half-year":
        synthesize_half_year(args.period)


if __name__ == "__main__":
    main()
