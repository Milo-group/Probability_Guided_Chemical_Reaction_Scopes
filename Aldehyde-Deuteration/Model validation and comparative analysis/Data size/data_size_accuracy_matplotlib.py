"""Plot data-size evaluation summary CSVs produced by Progressive_DataSize_Evaluation.R."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent

METRIC_COLS = [
    "train_accuracy",
    "test_accuracy",
    "external_accuracy",
    "smallest_group_recall",
]

COLORS = [
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#66A61E",
    "#E7298A",
    "#A6761D",
]

METRIC_COLORS = {col: COLORS[i % len(COLORS)] for i, col in enumerate(METRIC_COLS)}

SERIES_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]
SERIES_LINESTYLES: list[str | tuple[float, tuple[float, ...]]] = [
    "-",
    "--",
    "-.",
    ":",
    (0, (3, 1, 1, 1)),
    (0, (5, 2, 1, 2)),
    (0, (3, 3, 2, 2)),
]

FONT_TITLE = 16
FONT_AXIS = 14
FONT_TICK = 13
FONT_LEGEND = 13
FONT_LEGEND_PANEL = 12
MARKER_SIZE = 9
MARKER_SIZE_LEGEND = 11
LINE_WIDTH = 2.2
LINE_WIDTH_LEGEND = 2.8
MARKER_EDGE_WIDTH = 0.9
ERRORBAR_CAPSIZE = 4.5
ERRORBAR_LINEWIDTH = 1.4
FIG_HEIGHT_SIDE_BY_SIDE = 5.75
FIG_HEIGHT_OVERLAY = 5.5
FIG_WIDTH_PER_PANEL = 5.5
FIG_WIDTH_OVERLAY = 8.5

METRIC_LABELS = {
    "train_accuracy": "Train accuracy",
    "test_accuracy": "Internal validation accuracy",
    "external_accuracy": "External + unsampled accuracy",
    "smallest_group_recall": "Smallest-group recall",
}

METRIC_CLI_ALIASES = {
    "internal_validation_accuracy": "test_accuracy",
    "external_validation_accuracy": "external_accuracy",
    "train_accuracy": "train_accuracy",
    "smallest_group_recall": "smallest_group_recall",
    "test_accuracy": "test_accuracy",
}


def metric_color(metric_col: str) -> str:
    return METRIC_COLORS.get(metric_col, COLORS[0])


def series_marker_linestyle(series_index: int) -> tuple[str, str | tuple]:
    m = SERIES_MARKERS[series_index % len(SERIES_MARKERS)]
    ls = SERIES_LINESTYLES[series_index % len(SERIES_LINESTYLES)]
    return m, ls


def metric_label(metric_col: str) -> str:
    return METRIC_LABELS.get(metric_col, metric_col.replace("_", " "))


def resolve_metric_col(name: str) -> str:
    key = name.strip().lower().replace("-", "_")
    return METRIC_CLI_ALIASES.get(key, name)


def infer_raw_from_summary(summary_path: Path) -> Path:
    stem = summary_path.stem
    for old, new in (
        ("metrics_summary", "metrics_raw"),
        ("progressive_summary", "progressive_raw"),
        ("_summary", "_raw"),
    ):
        if old in stem:
            candidate = summary_path.with_name(stem.replace(old, new, 1) + summary_path.suffix)
            if candidate.is_file():
                return candidate
    return summary_path.with_name(f"{stem}_raw{summary_path.suffix}")


def resolve_path(path: Path, base_dir: Path) -> Path:
    return path if path.is_absolute() else base_dir / path


def build_panels(
    summaries: list[Path],
    raws: list[Path] | None,
    labels: list[str] | None,
    base_dir: Path,
) -> list[tuple[str, Path, Path]]:
    panels: list[tuple[str, Path, Path]] = []
    for i, summary in enumerate(summaries):
        sum_path = resolve_path(summary, base_dir)
        if raws and i < len(raws):
            raw_path = resolve_path(raws[i], base_dir)
        else:
            raw_path = infer_raw_from_summary(sum_path)
        if labels and i < len(labels):
            label = labels[i]
        else:
            label = f"Series {i + 1}"
        panels.append((label, sum_path, raw_path))
    return panels


def load_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)


def compute_z_based_error(
    summary: pd.DataFrame, col: str, repeat_counts: pd.Series
) -> pd.Series | None:
    sd_col = f"{col}_sd"
    if sd_col not in summary.columns:
        return None
    counts = summary["data_size"].map(repeat_counts).astype(float)
    z_95 = 1.96
    return z_95 * summary[sd_col].astype(float) / counts.apply(
        lambda n: math.sqrt(n) if pd.notna(n) and n > 0 else float("nan")
    )


def metric_legend_proxy_handles(summary: pd.DataFrame) -> list[Line2D]:
    handles: list[Line2D] = []
    for col in METRIC_COLS:
        if col not in summary.columns:
            continue
        c = metric_color(col)
        handles.append(
            Line2D(
                [0],
                [0],
                color=c,
                linestyle="-",
                marker="o",
                linewidth=LINE_WIDTH_LEGEND,
                markersize=MARKER_SIZE_LEGEND,
                markerfacecolor=c,
                markeredgecolor=c,
                label=metric_label(col),
            )
        )
    return handles


def plot_panel(
    ax,
    title: str,
    summary: pd.DataFrame,
    raw: pd.DataFrame,
    *,
    series_index: int = 0,
    show_legend: bool = True,
) -> None:
    x = summary["data_size"].values
    repeat_counts = raw.groupby("data_size")["repeat_id"].nunique()
    marker, linestyle = series_marker_linestyle(series_index)
    for col in METRIC_COLS:
        if col not in summary.columns:
            continue
        color = metric_color(col)
        yerr_series = compute_z_based_error(summary, col, repeat_counts)
        yerr = yerr_series.values if yerr_series is not None else None
        kw: dict = dict(
            color=color,
            marker=marker,
            linestyle=linestyle,
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=MARKER_EDGE_WIDTH,
            zorder=3,
        )
        if show_legend:
            kw["label"] = metric_label(col)
        ax.errorbar(
            x,
            summary[col].values,
            yerr=yerr,
            capsize=ERRORBAR_CAPSIZE,
            elinewidth=ERRORBAR_LINEWIDTH,
            capthick=ERRORBAR_LINEWIDTH,
            **kw,
        )

    ax.set_title(title, fontsize=FONT_TITLE, fontweight="semibold", pad=8)
    ax.set_xlabel("Sampled data size", fontsize=FONT_AXIS)
    ax.set_ylabel("Accuracy / recall", fontsize=FONT_AXIS)
    ax.set_ylim(0.5, 1.03)
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.tick_params(axis="both", labelsize=FONT_TICK)
    ax.grid(True, linestyle="-", alpha=0.25, linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_legend:
        ax.legend(loc="best", fontsize=FONT_LEGEND_PANEL, frameon=False, markerscale=1.2)


def save_side_by_side(
    panels: list[tuple[str, Path, Path]],
    out_path: Path,
) -> None:
    n = len(panels)
    fig_w = min(18.0, FIG_WIDTH_PER_PANEL * n)
    fig, axes = plt.subplots(1, n, figsize=(fig_w, FIG_HEIGHT_SIDE_BY_SIDE), sharey=True, squeeze=False)
    axes = axes.flatten()
    fig.subplots_adjust(left=0.085, right=0.98, bottom=0.42 if n >= 3 else 0.46, top=0.93, wspace=0.12)

    for si, (ax, (label, sum_path, raw_path)) in enumerate(zip(axes, panels)):
        plot_panel(
            ax,
            label,
            load_csv(sum_path),
            load_csv(raw_path),
            series_index=si,
            show_legend=False,
        )

    first_summary = load_csv(panels[0][1])
    proxy_handles = metric_legend_proxy_handles(first_summary)
    extra_artists: list = []
    if proxy_handles:
        legend = fig.legend(
            handles=proxy_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.30 if n >= 3 else 0.33),
            bbox_transform=fig.transFigure,
            ncol=min(4, len(proxy_handles)),
            fontsize=FONT_LEGEND,
            frameon=False,
            labelspacing=0.55,
            handlelength=2.2,
            handletextpad=0.65,
            columnspacing=1.1,
            borderaxespad=0.2,
            markerscale=1.15,
        )
        extra_artists.append(legend)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches="tight",
        bbox_extra_artists=extra_artists,
        pad_inches=0.18,
    )
    plt.close(fig)
    print(f"Saved: {out_path}")


def save_metric_overlay(
    panels: list[tuple[str, Path, Path]],
    metric_col: str,
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_OVERLAY, FIG_HEIGHT_OVERLAY))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.14, top=0.94)

    color = metric_color(metric_col)
    plotted = 0
    for si, (label, sum_path, raw_path) in enumerate(panels):
        summary = load_csv(sum_path)
        raw = load_csv(raw_path)
        if metric_col not in summary.columns:
            print(f"Skipping {label}: column {metric_col!r} not found", file=sys.stderr)
            continue
        marker, linestyle = series_marker_linestyle(si)
        repeat_counts = raw.groupby("data_size")["repeat_id"].nunique()
        yerr_series = compute_z_based_error(summary, metric_col, repeat_counts)
        yerr = yerr_series.values if yerr_series is not None else None
        ax.errorbar(
            summary["data_size"].values,
            summary[metric_col].values,
            yerr=yerr,
            capsize=ERRORBAR_CAPSIZE,
            elinewidth=ERRORBAR_LINEWIDTH,
            capthick=ERRORBAR_LINEWIDTH,
            color=color,
            marker=marker,
            linestyle=linestyle,
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=MARKER_EDGE_WIDTH,
            label=label,
            zorder=3,
        )
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        raise KeyError(f"No series contained column {metric_col!r}")

    ax.set_xlabel("Sampled data size", fontsize=FONT_AXIS)
    ax.set_ylabel(metric_label(metric_col), fontsize=FONT_AXIS)
    ax.set_ylim(0.5, 1.03)
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.tick_params(axis="both", labelsize=FONT_TICK)
    ax.grid(True, linestyle="-", alpha=0.25, linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", fontsize=FONT_LEGEND, frameon=False, markerscale=1.2)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print(f"Saved: {out_path}")


def slugify(text: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "_", text.lower()).strip("_")
    return s[:24] if s else "series"


def run_plots(
    panels: list[tuple[str, Path, Path]],
    *,
    metrics: list[str],
    out_dir: Path,
    out_prefix: str,
) -> None:
    if not panels:
        raise SystemExit("Provide at least one --summary file")

    slug = "_vs_".join(slugify(label) for label, _, _ in panels[:4])
    if len(panels) > 4:
        slug = f"{slug}_and_{len(panels) - 4}_more"

    save_side_by_side(panels, out_dir / f"{out_prefix}_{slug}.png")

    if len(panels) < 2:
        return

    for name in metrics:
        metric_col = resolve_metric_col(name)
        if metric_col not in METRIC_COLS:
            print(f"Skipping unknown metric: {name!r}", file=sys.stderr)
            continue
        try:
            save_metric_overlay(
                panels,
                metric_col,
                out_dir / f"{out_prefix}_{slug}_{metric_col}.png",
            )
        except (KeyError, OSError) as exc:
            print(exc, file=sys.stderr)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot data-size evaluation summary tables.")
    p.add_argument(
        "--summary",
        nargs="+",
        type=Path,
        required=True,
        help="Summary table from the R evaluation (one per condition to compare)",
    )
    p.add_argument(
        "--raw",
        nargs="*",
        type=Path,
        help="Matching raw tables (optional; paired raw is inferred when omitted)",
    )
    p.add_argument("--labels", nargs="+", help="Legend / panel titles (same order as --summary)")
    p.add_argument(
        "--metrics",
        nargs="+",
        default=["external_validation_accuracy", "smallest_group_recall"],
        help="Metrics for per-metric overlay figures when comparing 2+ series",
    )
    p.add_argument("--out-prefix", default="data_size_plot", help="Output figure name prefix")
    p.add_argument(
        "--data-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Base directory for relative --summary / --raw paths",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for output figures (default: same as --data-dir)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    out_dir = (args.output_dir or args.data_dir).resolve()

    if args.raw and len(args.raw) not in (0, len(args.summary)):
        raise SystemExit("--raw must be omitted or provide one path per --summary")

    panels = build_panels(args.summary, args.raw or None, args.labels, data_dir)
    run_plots(
        panels,
        metrics=args.metrics,
        out_dir=out_dir,
        out_prefix=args.out_prefix,
    )


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
