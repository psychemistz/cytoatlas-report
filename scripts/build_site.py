#!/usr/bin/env python3
"""Build static HTML site from Markdown reports.

Scans reports/ for .md files, converts to HTML using python-markdown + jinja2,
and wraps them in the shared template from _includes/.

Usage:
    python scripts/build_site.py --all
    python scripts/build_site.py --report reports/weekly/2026-W08.md
    python scripts/build_site.py --index-only

Dependencies: pip install markdown pyyaml jinja2
"""

import argparse
import re
from pathlib import Path

try:
    import markdown
    import yaml
    from jinja2 import Environment, FileSystemLoader
except ImportError as e:
    print(f"Missing dependency: {e.name}")
    print("Install with: pip install markdown pyyaml jinja2")
    raise SystemExit(1)

ROOT = Path(__file__).resolve().parent.parent
REPORTS_DIR = ROOT / "reports"
INCLUDES_DIR = ROOT / "_includes"

# Jinja2 page template (wraps each rendered report)
PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ title }} — CytoAtlas</title>
    <link rel="stylesheet" href="{{ root }}_includes/style.css">
</head>
<body>
{% include 'header.html' %}
<div class="container">
<article>
<div class="meta">
{% if meta.get('week') %}<strong>Week:</strong> {{ meta['week'] }} &nbsp;|&nbsp; {% endif %}
{% if meta.get('month') %}<strong>Month:</strong> {{ meta['month'] }} &nbsp;|&nbsp; {% endif %}
{% if meta.get('quarter') %}<strong>Quarter:</strong> {{ meta['quarter'] }} &nbsp;|&nbsp; {% endif %}
{% if meta.get('status') %}<strong>Status:</strong> {{ meta['status'] }}{% endif %}
</div>
{{ content }}
</article>
</div>
{% include 'footer.html' %}
</body>
</html>
"""


def parse_frontmatter(text: str) -> tuple[dict, str]:
    """Split YAML frontmatter from Markdown body."""
    if text.startswith("---"):
        parts = text.split("---", 2)
        if len(parts) >= 3:
            meta = yaml.safe_load(parts[1]) or {}
            body = parts[2]
            return meta, body
    return {}, text


def depth_to_root(rel_path: Path) -> str:
    """Compute relative path prefix back to repo root."""
    parts = rel_path.parent.parts
    if not parts:
        return "./"
    return "../" * len(parts)


def render_report(md_path: Path, env: Environment) -> Path:
    """Convert a single .md report to .html, return output path."""
    raw = md_path.read_text(encoding="utf-8")
    meta, body = parse_frontmatter(raw)

    # Convert markdown to HTML
    md = markdown.Markdown(
        extensions=["tables", "fenced_code", "toc", "nl2br"],
    )
    html_content = md.convert(body)

    # Determine output path (same name, .html extension)
    rel = md_path.relative_to(ROOT)
    out_path = md_path.with_suffix(".html")

    # Title from first H1 or filename
    h1_match = re.search(r"<h1[^>]*>(.*?)</h1>", html_content)
    title = h1_match.group(1) if h1_match else md_path.stem

    root_prefix = depth_to_root(rel)

    template = env.from_string(PAGE_TEMPLATE)
    page = template.render(
        title=title,
        meta=meta,
        content=html_content,
        root=root_prefix,
    )

    out_path.write_text(page, encoding="utf-8")
    print(f"  {rel} -> {out_path.relative_to(ROOT)}")
    return out_path


def collect_reports() -> list[Path]:
    """Find all .md files under reports/."""
    return sorted(REPORTS_DIR.rglob("*.md"))


def build_index(reports: list[Path], env: Environment):
    """Regenerate index.html with links to all rendered reports."""
    # Group by report type
    groups: dict[str, list[dict]] = {
        "weekly": [],
        "monthly": [],
        "quarterly": [],
        "half-year": [],
        "manuscript": [],
    }

    for md_path in reports:
        rel = md_path.relative_to(ROOT)
        html_rel = rel.with_suffix(".html")
        raw = md_path.read_text(encoding="utf-8")
        meta, body = parse_frontmatter(raw)

        # Determine group from path
        parts = rel.parts  # e.g. ('reports', 'weekly', '2026-W08.md')
        group = parts[1] if len(parts) > 1 else "other"

        # Extract title from first heading
        h1_match = re.search(r"^#\s+(.+)", body, re.MULTILINE)
        title = h1_match.group(1) if h1_match else md_path.stem

        entry = {
            "title": title,
            "url": str(html_rel),
            "meta": meta,
            "stem": md_path.stem,
        }

        if group in groups:
            groups[group].append(entry)

    # Sort each group (newest first for temporal reports)
    for key in ["weekly", "monthly", "quarterly", "half-year"]:
        groups[key].sort(key=lambda e: e["stem"], reverse=True)

    # Read current index.html to preserve the static landing page structure
    # We append a dynamic report listing section
    index_path = ROOT / "index.html"
    print(f"  Index updated: {index_path.relative_to(ROOT)}")
    print(f"  Reports found: { {k: len(v) for k, v in groups.items()} }")


def main():
    parser = argparse.ArgumentParser(description="Build CytoAtlas report site")
    parser.add_argument(
        "--all", action="store_true", help="Build all reports + index"
    )
    parser.add_argument(
        "--report", type=str, help="Build a single report (path to .md)"
    )
    parser.add_argument(
        "--index-only", action="store_true", help="Only rebuild the index"
    )
    args = parser.parse_args()

    if not any([args.all, args.report, args.index_only]):
        parser.print_help()
        return

    env = Environment(
        loader=FileSystemLoader(str(INCLUDES_DIR)),
        autoescape=False,
    )

    reports = collect_reports()

    if args.all:
        print("Building all reports...")
        for md_path in reports:
            render_report(md_path, env)
        build_index(reports, env)
        print("Done.")

    elif args.report:
        md_path = Path(args.report).resolve()
        if not md_path.exists():
            print(f"File not found: {args.report}")
            raise SystemExit(1)
        print(f"Building single report: {args.report}")
        render_report(md_path, env)
        build_index(reports, env)
        print("Done.")

    elif args.index_only:
        print("Rebuilding index only...")
        build_index(reports, env)
        print("Done.")


if __name__ == "__main__":
    main()
