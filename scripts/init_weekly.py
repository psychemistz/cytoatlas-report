#!/usr/bin/env python3
"""Scaffold a new weekly report from the template.

Usage:
    python scripts/init_weekly.py              # current ISO week
    python scripts/init_weekly.py --week 2026-W08
"""

import argparse
import re
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "templates" / "weekly.md"
REPORTS_DIR = ROOT / "reports" / "weekly"
FIGURES_DIR = ROOT / "figures" / "weekly"


def iso_week_str(dt: datetime) -> str:
    year, week, _ = dt.isocalendar()
    return f"{year}-W{week:02d}"


def week_date_range(week_str: str) -> tuple[str, str]:
    """Return (monday, sunday) date strings for a given YYYY-WXX."""
    year, week = int(week_str[:4]), int(week_str.split("W")[1])
    # ISO: week 1 contains the year's first Thursday
    jan4 = datetime(year, 1, 4)
    start_of_w1 = jan4 - timedelta(days=jan4.weekday())
    monday = start_of_w1 + timedelta(weeks=week - 1)
    sunday = monday + timedelta(days=6)
    return monday.strftime("%Y-%m-%d"), sunday.strftime("%Y-%m-%d")


def previous_week_next_steps(week_str: str) -> str | None:
    """Extract 'Next Steps' from the previous week's report, if it exists."""
    year, week = int(week_str[:4]), int(week_str.split("W")[1])
    if week == 1:
        prev = f"{year - 1}-W52"
    else:
        prev = f"{year}-W{week - 1:02d}"

    prev_file = REPORTS_DIR / f"{prev}.md"
    if not prev_file.exists():
        return None

    text = prev_file.read_text()
    match = re.search(
        r"## Next Steps\s*\n(.*?)(?=\n## |\Z)", text, re.DOTALL
    )
    if match:
        return match.group(1).strip()
    return None


def main():
    parser = argparse.ArgumentParser(description="Initialize a weekly report")
    parser.add_argument(
        "--week",
        type=str,
        default=None,
        help="Week in YYYY-WXX format (default: current week)",
    )
    args = parser.parse_args()

    week = args.week or iso_week_str(datetime.now())

    # Validate format
    if not re.match(r"^\d{4}-W\d{2}$", week):
        parser.error(f"Invalid week format: {week!r}. Use YYYY-WXX.")

    report_path = REPORTS_DIR / f"{week}.md"
    figures_path = FIGURES_DIR / week

    if report_path.exists():
        print(f"Report already exists: {report_path}")
        return

    # Read template
    template = TEMPLATE.read_text()

    # Fill frontmatter
    mon, sun = week_date_range(week)
    content = template.replace("YYYY-WXX", week)
    content = content.replace("YYYY-MM-DD to YYYY-MM-DD", f"{mon} to {sun}")
    content = content.replace('week: "{{ week }}"', f'week: "{week}"')

    # Replace Jinja-style placeholders in body
    content = content.replace("{{ week }}", week)
    content = content.replace("{{ date_range }}", f"{mon} to {sun}")

    # Create report and figures directory
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    figures_path.mkdir(parents=True, exist_ok=True)

    report_path.write_text(content)
    print(f"Created: {report_path}")
    print(f"Created: {figures_path}/")

    # Show previous week's next steps as reminder
    prev_steps = previous_week_next_steps(week)
    if prev_steps:
        print(f"\n--- Last week's Next Steps ---")
        print(prev_steps)
        print("---")
    else:
        print("\nNo previous weekly report found.")


if __name__ == "__main__":
    main()
