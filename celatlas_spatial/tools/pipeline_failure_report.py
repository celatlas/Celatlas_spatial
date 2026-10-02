#!/usr/bin/env python3
"""Generate a minimal HTML report when the Celatlas pipeline fails.

The regular report generators need enough downstream outputs to render a full
analysis report. This fallback uses only the sample directory and pipeline log,
so production callers still get a stable HTML artifact after early failures.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import html
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def _read_tail(path: Path, max_lines: int) -> List[str]:
    if not path or not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()
    except OSError as exc:
        return [f"Failed to read log file {path}: {exc}"]
    return [line.rstrip("\n") for line in lines[-max_lines:]]


def _exists_text(path: Path) -> str:
    return "available" if path.exists() else "missing"


def _stage_candidates(sample_dir: Path, sample: str, mode: str) -> List[Tuple[str, Path]]:
    common = [
        ("00.sample", sample_dir / "00.sample"),
        ("01.barcode", sample_dir / "01.barcode"),
        ("02.cutadapt", sample_dir / "02.cutadapt"),
        ("03.star", sample_dir / "03.star"),
        ("04.featureCounts", sample_dir / "04.featureCounts"),
        ("05.count", sample_dir / "05.count"),
    ]
    if mode == "scrna":
        common.extend(
            [
                ("05.count matrix", sample_dir / "05.count" / f"{sample}_filtered_feature_bc_matrix"),
                ("06_analysis_wrapper", sample_dir / "06_analysis_wrapper"),
            ]
        )
    else:
        common.extend(
            [
                ("06.segment mask", sample_dir / "06.segment" / "mask"),
                ("06.segment binSegment", sample_dir / "06.segment" / "01.binsegment"),
                ("06.segment square_bin", sample_dir / "06.segment" / "01.binsegment" / "square_bin"),
                ("07.outs binned outputs", sample_dir / "07.outs" / "binned_outputs"),
            ]
        )
    return common


def _list_artifacts(sample_dir: Path, sample: str, mode: str) -> List[Path]:
    names = [
        sample_dir / "pipeline.log",
        sample_dir / "reanalysis_pipeline.log",
        sample_dir / ".fastq_inputs.tsv",
        sample_dir / "01.barcode" / f"{sample}_summary.txt",
        sample_dir / "01.barcode" / f"{sample}_barcode_stat.txt",
        sample_dir / "05.count" / f"{sample}_count_detail.txt",
    ]
    if mode == "scrna":
        names.extend(
            [
                sample_dir / "05.count" / f"{sample}_filtered_feature_bc_matrix",
                sample_dir / "06_analysis_wrapper",
            ]
        )
    else:
        names.extend(
            [
                sample_dir / "06.segment" / "01.binsegment" / "images",
                sample_dir / "06.segment" / "01.binsegment" / "square_bin",
                sample_dir / "07.outs" / "binned_outputs",
            ]
        )
    return [path for path in names if path.exists()]


def _rel(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def _write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def _html_list(items: Iterable[str]) -> str:
    return "\n".join(f"<li>{item}</li>" for item in items)


def _normalize_metric_display(value: str) -> str:
    value = value.strip()
    value = re.sub(r"([0-9,])\(", r"\1 (", value)
    return value


def _format_count(value: str) -> str:
    raw = str(value).replace(",", "").strip()
    if re.fullmatch(r"-?\d+", raw):
        return f"{int(raw):,}"
    return str(value).strip()


def _parse_step_log(path: Path) -> Dict[str, Dict[str, str]]:
    result = {"arguments": {}, "metrics": {}}
    if not path.exists():
        return result

    section = ""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                stripped = line.strip()
                if stripped in {"Arguments", "Metrics", "Hidden Metrics"}:
                    section = stripped
                    continue
                if not stripped or set(stripped) == {"-"}:
                    continue
                if stripped in {
                    "Commands",
                    "Help Content",
                    "Data Summary",
                    "Summary Files",
                    "Output Directory Entries",
                }:
                    section = ""
                    continue

                if section == "Arguments" and ":" in line and not line.startswith(" "):
                    key, value = line.split(":", 1)
                    result["arguments"][key.strip()] = value.strip()
                elif section in {"Metrics", "Hidden Metrics"}:
                    if not line.startswith("- ") or ":" not in line:
                        if not line.startswith("  "):
                            section = ""
                        continue
                    key, value = line[2:].split(":", 1)
                    result["metrics"][key.strip()] = _normalize_metric_display(value)
    except OSError as exc:
        result["metrics"]["Step Log Parse Error"] = str(exc)
    return result


def _parse_stat_file(path: Path) -> Dict[str, str]:
    stats: Dict[str, str] = {}
    if not path.exists():
        return stats
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if ":" not in line:
                    continue
                key, value = line.split(":", 1)
                stats[key.strip()] = _normalize_metric_display(value)
    except OSError:
        return stats
    return stats


def _parse_star_final_log(path: Path) -> Dict[str, str]:
    stats: Dict[str, str] = {}
    if not path.exists():
        return stats

    key_map = {
        "Number of input reads": "STAR Input Reads",
        "Uniquely mapped reads number": "Uniquely Mapped Reads",
        "Uniquely mapped reads %": "Uniquely Mapped Reads Percentage",
        "Number of reads mapped to multiple loci": "Multi-Mapped Reads",
        "% of reads mapped to multiple loci": "Multi-Mapped Reads Percentage",
        "Number of reads mapped to too many loci": "Reads Mapped To Too Many Loci",
        "% of reads mapped to too many loci": "Reads Mapped To Too Many Loci Percentage",
        "Number of reads unmapped: too short": "Unmapped Reads - Too Short",
        "% of reads unmapped: too short": "Unmapped Reads - Too Short Percentage",
        "Number of reads unmapped: other": "Unmapped Reads - Other",
        "% of reads unmapped: other": "Unmapped Reads - Other Percentage",
    }
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "|" not in line:
                    continue
                key, value = [item.strip() for item in line.split("|", 1)]
                mapped = key_map.get(key)
                if not mapped:
                    continue
                value = value.replace("%", "").strip() if mapped.endswith("Percentage") else value.strip()
                stats[mapped] = f"{value}%" if mapped.endswith("Percentage") else _format_count(value)
    except OSError:
        return stats
    return stats


def _metric(metrics: Dict[str, str], key: str, *aliases: str) -> str:
    for candidate in (key, *aliases):
        value = metrics.get(candidate)
        if value:
            if "(" not in value:
                pct = metrics.get(f"{candidate} Percentage")
                if pct:
                    return f"{_format_count(value)} ({pct})"
            return _normalize_metric_display(value)
    return ""


def _read_downsample(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except OSError:
        return []


def _bin_sort_key(label: str) -> Tuple[int, str]:
    normalized = label.lower().replace("bin", "")
    if normalized == "raw":
        return (10_000, label)
    if normalized.isdigit():
        return (int(normalized), label)
    return (9_000, label)


def _bin_label(path: Path, sample: str) -> str:
    name = path.parent.name
    prefix = f"{sample}_"
    if name.startswith(prefix):
        name = name[len(prefix):]
    if name.lower().startswith("bin"):
        return name[3:]
    return name


def _collect_bin_metrics(sample_dir: Path, sample: str) -> List[Dict[str, str]]:
    square_bin_roots = [
        sample_dir / "06.segment" / "01.binsegment" / "square_bin",
    ]
    rows: List[Dict[str, str]] = []
    seen = set()
    wanted = [
        "Total Genes",
        "Estimated Number of square bin",
        "Fraction Reads in Square",
        "Mean Reads per square bin",
        "Median Reads per square bin",
        "Mean UMI per square bin",
        "Median UMI per square bin",
        "Mean Genes per square bin",
        "Median Genes per square bin",
        "Cavity Filter Removed Fraction",
    ]
    for root in square_bin_roots:
        if not root.exists():
            continue
        for stat_file in sorted(root.glob(f"{sample}_*/stat.txt")):
            label = _bin_label(stat_file, sample)
            if label in seen:
                continue
            seen.add(label)
            stats = _parse_stat_file(stat_file)
            if not stats:
                continue
            row = {"Bin Size": label}
            for key in wanted:
                row[key] = stats.get(key, "")
            rows.append(row)
    return sorted(rows, key=lambda item: _bin_sort_key(item["Bin Size"]))


def _collect_qc_summary(sample_dir: Path, sample: str) -> Dict:
    barcode = _parse_step_log(sample_dir / "01.barcode" / "step.log")
    star = _parse_step_log(sample_dir / "03.star" / "step.log")
    star_final = _parse_star_final_log(sample_dir / "03.star" / f"{sample}_Log.final.out")
    feature = _parse_step_log(sample_dir / "04.featureCounts" / "step.log")
    count = _parse_step_log(sample_dir / "05.count" / "step.log")

    star_metrics = {**star_final, **star["metrics"]}
    genome_dir = star["arguments"].get("genomeDir", "")
    genome = Path(genome_dir).name if genome_dir else _metric(star_metrics, "Genome")

    feature_type = _metric(feature["metrics"], "Feature Type")
    feature_type = feature_type.lower() if feature_type else ""

    return {
        "barcode": [
            ("Number of Reads", _metric(barcode["metrics"], "Raw Reads")),
            ("Valid Reads", _metric(barcode["metrics"], "Valid Reads")),
            ("Exactly Matched Reads", _metric(barcode["metrics"], "Exactly Matched Reads")),
            ("Mismatched Reads", _metric(barcode["metrics"], "Mismatched Reads")),
            ("Q30 of Barcodes", _metric(barcode["metrics"], "Q30 of Barcodes")),
            ("Q30 of UMIs", _metric(barcode["metrics"], "Q30 of UMIs")),
        ],
        "mapping": [
            ("Genome", genome),
            ("STAR Input Reads", _metric(star_metrics, "STAR Input Reads")),
            ("Uniquely Mapped Reads", _metric(star_metrics, "Uniquely Mapped Reads")),
            ("Multi-Mapped Reads", _metric(star_metrics, "Multi-Mapped Reads")),
            (
                "Mapped to Too Many Loci",
                _metric(
                    star_metrics,
                    "Reads Mapped To Too Many Loci",
                    "Reads Mapped to Too Many Loci",
                ),
            ),
            ("Unmapped Reads - Too Short", _metric(star_metrics, "Unmapped Reads - Too Short")),
            ("Unmapped Reads - Other", _metric(star_metrics, "Unmapped Reads - Other")),
        ],
        "featurecounts": [
            ("Feature Type", feature_type),
            ("Reads Assigned To Exonic Regions", _metric(feature["metrics"], "Reads Assigned To Exonic Regions")),
            ("Reads Assigned To Intronic Regions", _metric(feature["metrics"], "Reads Assigned To Intronic Regions")),
            ("Reads Assigned To Intergenic Regions", _metric(feature["metrics"], "Reads Assigned To Intergenic Regions")),
            ("Reads Unassigned Ambiguity", _metric(feature["metrics"], "Reads Unassigned Ambiguity")),
        ],
        "count": [
            ("Estimated Number of Cells", _metric(count["metrics"], "Estimated Number of Cells")),
            ("Fraction Reads in Cells", _metric(count["metrics"], "Fraction Reads in Cells")),
            ("Mean Reads per Cell", _metric(count["metrics"], "Mean Reads per Cell")),
            ("Median UMI per Cell", _metric(count["metrics"], "Median UMI per Cell")),
            ("Total Genes", _metric(count["metrics"], "Total Genes")),
            ("Median Genes per Cell", _metric(count["metrics"], "Median Genes per Cell")),
            ("Saturation", _metric(count["metrics"], "Saturation")),
        ],
        "downsample": _read_downsample(sample_dir / "05.count" / f"{sample}_downsample.tsv"),
        "bin_metrics": _collect_bin_metrics(sample_dir, sample),
    }


def _has_rows(rows: List[Tuple[str, str]]) -> bool:
    return any(value for _, value in rows)


def _render_kv_table(rows: List[Tuple[str, str]]) -> str:
    rendered = []
    for key, value in rows:
        if not value:
            continue
        rendered.append(
            "<tr>"
            f"<td>{html.escape(key)}</td>"
            f"<td>{html.escape(str(value))}</td>"
            "</tr>"
        )
    return "\n".join(rendered)


def _render_dict_table(rows: List[Dict[str, str]], columns: List[Tuple[str, str]]) -> str:
    if not rows:
        return ""
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key, _ in columns)
            + "</tr>"
        )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _render_qc_sections(qc_summary: Dict) -> str:
    sections = []
    table_specs = [
        ("Barcode QC", qc_summary.get("barcode", [])),
        ("Mapping", qc_summary.get("mapping", [])),
        ("FeatureCounts", qc_summary.get("featurecounts", [])),
        ("Count And Saturation Summary", qc_summary.get("count", [])),
    ]
    for title, rows in table_specs:
        if not _has_rows(rows):
            continue
        sections.append(
            '<section class="panel" style="margin-top:18px;">'
            f"<h2>{html.escape(title)}</h2>"
            "<table><tbody>"
            f"{_render_kv_table(rows)}"
            "</tbody></table>"
            "</section>"
        )

    downsample = qc_summary.get("downsample", [])
    if downsample:
        sections.append(
            '<section class="panel" style="margin-top:18px;">'
            "<h2>Sequencing Saturation</h2>"
            + _render_dict_table(
                downsample,
                [
                    ("read_fraction", "Read Fraction"),
                    ("median_gene_number", "Median Genes"),
                    ("umi_saturation", "UMI Saturation (%)"),
                    ("read_saturation", "Read Saturation (%)"),
                ],
            )
            + "</section>"
        )

    bin_metrics = qc_summary.get("bin_metrics", [])
    if bin_metrics:
        sections.append(
            '<section class="panel" style="margin-top:18px;">'
            "<h2>Bin Metrics Overview</h2>"
            + _render_dict_table(
                bin_metrics,
                [
                    ("Bin Size", "Bin Size"),
                    ("Total Genes", "Total Genes"),
                    ("Estimated Number of square bin", "Estimated Squares"),
                    ("Fraction Reads in Square", "Fraction Reads"),
                    ("Mean Reads per square bin", "Mean Reads"),
                    ("Median Reads per square bin", "Median Reads"),
                    ("Mean UMI per square bin", "Mean UMI"),
                    ("Median UMI per square bin", "Median UMI"),
                    ("Mean Genes per square bin", "Mean Genes"),
                    ("Median Genes per square bin", "Median Genes"),
                    ("Cavity Filter Removed Fraction", "Cavity Removed"),
                ],
            )
            + "</section>"
        )
    return "\n".join(sections)


def build_html(payload: Dict, log_tail: List[str], stages: List[Tuple[str, Path]], artifacts: List[Path]) -> str:
    sample_dir = Path(payload["sample_dir"])
    qc_html = _render_qc_sections(payload.get("qc_summary", {}))
    partial_report = payload.get("partial_report") or ""
    partial_html = ""
    if partial_report:
        partial_path = Path(partial_report)
        if partial_path.exists():
            partial_html = (
                '<div class="callout">'
                "<strong>Partial analysis report:</strong> "
                f'<a href="{html.escape(_rel(partial_path, sample_dir))}">{html.escape(partial_path.name)}</a>'
                "</div>"
            )

    stage_rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f"<td><span class=\"status {_exists_text(path)}\">{html.escape(_exists_text(path))}</span></td>"
        f"<td>{html.escape(_rel(path, sample_dir))}</td>"
        "</tr>"
        for name, path in stages
    )

    artifact_items = []
    for path in artifacts:
        label = html.escape(_rel(path, sample_dir))
        if path.is_file():
            artifact_items.append(f'<a href="{label}">{label}</a>')
        else:
            artifact_items.append(label)

    log_text = html.escape("\n".join(log_tail) if log_tail else "No pipeline log was available.")
    generated_at = html.escape(payload["generated_at"])
    failed_command = html.escape(payload.get("failed_command") or "unknown")
    workflow = html.escape(payload.get("workflow") or "unknown")
    report_kind = html.escape(payload.get("report_kind") or "unknown")

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(payload["sample"])} Pipeline Failure Report</title>
  <style>
    :root {{
      --ink: #1f2937;
      --muted: #667085;
      --line: #d0d5dd;
      --panel: #ffffff;
      --bg: #f6f8fb;
      --danger: #b42318;
      --danger-bg: #fff1f0;
      --ok: #067647;
      --ok-bg: #ecfdf3;
      --warn: #b54708;
      --warn-bg: #fffaeb;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: Arial, Helvetica, sans-serif;
      line-height: 1.45;
    }}
    .wrap {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 28px;
    }}
    .hero {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-left: 6px solid var(--danger);
      padding: 24px;
      margin-bottom: 18px;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 28px;
      letter-spacing: 0;
    }}
    h2 {{
      margin: 0 0 14px;
      font-size: 18px;
      letter-spacing: 0;
    }}
    .sub {{
      color: var(--muted);
      margin: 0;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin: 18px 0;
    }}
    .metric, .panel, .callout {{
      background: var(--panel);
      border: 1px solid var(--line);
      padding: 16px;
    }}
    .metric .label {{
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
    }}
    .metric .value {{
      margin-top: 6px;
      font-size: 18px;
      font-weight: 700;
      overflow-wrap: anywhere;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: var(--panel);
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      padding: 10px;
      text-align: left;
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
    }}
    .status {{
      display: inline-block;
      min-width: 74px;
      text-align: center;
      padding: 3px 8px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
    }}
    .status.available {{
      background: var(--ok-bg);
      color: var(--ok);
    }}
    .status.missing {{
      background: var(--warn-bg);
      color: var(--warn);
    }}
    .callout {{
      background: var(--danger-bg);
      border-color: #fecdca;
      color: var(--danger);
      margin-bottom: 18px;
    }}
    pre {{
      margin: 0;
      padding: 14px;
      overflow: auto;
      background: #111827;
      color: #f9fafb;
      font-size: 12px;
      max-height: 520px;
      white-space: pre-wrap;
    }}
    ul {{
      margin-top: 0;
    }}
    a {{
      color: #175cd3;
    }}
  </style>
</head>
<body>
  <main class="wrap">
    <section class="hero">
      <h1>Celatlas Pipeline Failure Report</h1>
      <p class="sub">The pipeline stopped before a complete {report_kind} analysis report could be produced. This diagnostic report was generated automatically.</p>
    </section>

    <div class="callout">
      <strong>Failure command:</strong> <code>{failed_command}</code>
    </div>
    {partial_html}

    <section class="grid">
      <div class="metric"><div class="label">Sample</div><div class="value">{html.escape(payload["sample"])}</div></div>
      <div class="metric"><div class="label">Workflow</div><div class="value">{workflow}</div></div>
      <div class="metric"><div class="label">Mode</div><div class="value">{html.escape(payload.get("mode") or "unknown")}</div></div>
      <div class="metric"><div class="label">Method</div><div class="value">{html.escape(payload.get("method") or "unknown")}</div></div>
      <div class="metric"><div class="label">Chemistry</div><div class="value">{html.escape(payload.get("chemistry") or "unknown")}</div></div>
      <div class="metric"><div class="label">Species</div><div class="value">{html.escape(payload.get("species") or "unknown")}</div></div>
      <div class="metric"><div class="label">Exit Code</div><div class="value">{html.escape(str(payload.get("exit_code", "unknown")))}</div></div>
      <div class="metric"><div class="label">Generated</div><div class="value">{generated_at}</div></div>
    </section>

    {qc_html}

    <section class="panel">
      <h2>Stage Availability</h2>
      <table>
        <thead><tr><th>Stage</th><th>Status</th><th>Path</th></tr></thead>
        <tbody>{stage_rows}</tbody>
      </table>
    </section>

    <section class="panel" style="margin-top:18px;">
      <h2>Available Artifacts</h2>
      <ul>{_html_list(artifact_items) if artifact_items else "<li>No standard artifacts were detected yet.</li>"}</ul>
    </section>

    <section class="panel" style="margin-top:18px;">
      <h2>Pipeline Log Tail</h2>
      <pre>{log_text}</pre>
    </section>
  </main>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a Celatlas pipeline failure report")
    parser.add_argument("--sampledir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", default="")
    parser.add_argument("--method", default="")
    parser.add_argument("--chemistry", default="")
    parser.add_argument("--species", default="")
    parser.add_argument("--workflow", default="")
    parser.add_argument("--log-file", default="")
    parser.add_argument("--exit-code", default="")
    parser.add_argument("--line-no", default="")
    parser.add_argument("--failed-command", default="")
    parser.add_argument("--report-kind", choices=["scrna", "spatial"], default="spatial")
    parser.add_argument("--output-filename", required=True)
    parser.add_argument("--partial-report", default="")
    parser.add_argument("--log-lines", type=int, default=220)
    args = parser.parse_args()

    sample_dir = Path(args.sampledir).resolve()
    sample_dir.mkdir(parents=True, exist_ok=True)
    output_path = sample_dir / args.output_filename
    failure_report_path = sample_dir / f"{args.sample}_pipeline_failure_report.html"
    summary_path = sample_dir / f"{args.sample}_pipeline_failure_summary.json"
    log_file = Path(args.log_file).resolve() if args.log_file else Path()

    payload = {
        "sample_dir": str(sample_dir),
        "sample": args.sample,
        "mode": args.mode,
        "method": args.method,
        "chemistry": args.chemistry,
        "species": args.species,
        "workflow": args.workflow,
        "log_file": str(log_file) if log_file else "",
        "exit_code": args.exit_code,
        "line_no": args.line_no,
        "failed_command": args.failed_command,
        "report_kind": args.report_kind,
        "output_report": str(output_path),
        "failure_report": str(failure_report_path),
        "partial_report": args.partial_report,
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
    }

    stages = _stage_candidates(sample_dir, args.sample, args.mode)
    artifacts = _list_artifacts(sample_dir, args.sample, args.mode)
    payload["stages"] = [
        {"name": name, "path": str(path), "status": _exists_text(path)}
        for name, path in stages
    ]
    payload["artifacts"] = [str(path) for path in artifacts]
    payload["qc_summary"] = _collect_qc_summary(sample_dir, args.sample)

    log_tail = _read_tail(log_file, args.log_lines)
    payload["log_tail"] = log_tail

    html_text = build_html(payload, log_tail, stages, artifacts)
    for path in (output_path, failure_report_path):
        with path.open("w", encoding="utf-8") as handle:
            handle.write(html_text)
    _write_json(summary_path, payload)

    print(f"Failure report generated: {output_path}")
    print(f"Failure summary generated: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
