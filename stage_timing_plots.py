"""Pie charts for paired FGO phase timings, without rerunning estimators."""

import numpy as np


FIELDS = ("add_state_factor_mean_ms", "estimate_mean_ms", "marginalize_mean_ms")
COLORS = ("#159C99", "#2F6BA5", "#DF9748")
METHODS = {"schur": "Schur complement", "discard": "Direct discard"}


def values(row):
    durations = np.array([row[field] for field in FIELDS])
    return durations, durations / durations.sum() * 100.0


def percent(value, digits=1):
    threshold = 10 ** -digits
    return f"<{threshold:.{digits}f}%" if 0 < value < threshold else f"{value:.{digits}f}%"


def render_stage_plots(output, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.fonttype": "none"})
    lookup = {(r["window"], r["policy"]): r for r in summary["aggregates"]}
    windows = summary["windows"]
    settings = (f"Seed 7 · max {summary['max_iteration']} GN iterations/update · BLAS: 1 thread · "
                f"outlier weight {summary['data']['outlier_weight']:g}")
    behavior = ("Local Schur prior; W=100 retains the complete graph. Detailed timings are available in CSV."
                if summary.get("marginalization_algorithm") else
                "Existing algorithm retained, including initial-state removal at W=100. Detailed timings are available in CSV.")

    def pie(ax, row, threshold=5, fontsize=14):
        durations, percentages = values(row)
        ax.pie(durations, colors=COLORS, startangle=90, counterclock=False,
               autopct=lambda pct: f"{pct:.1f}%" if pct >= threshold else "",
               pctdistance=.66, radius=1,
               wedgeprops={"edgecolor": "white", "linewidth": 1.2},
               textprops={"fontsize": fontsize, "fontweight": "bold", "color": "white"})
        ax.set_aspect("equal")
        return durations, percentages

    def save(fig, name):
        fig.savefig(output / f"{name}.png", dpi=180, facecolor="white")
        fig.savefig(output / f"{name}.svg", facecolor="white")
        plt.close(fig)

    for window in windows:
        fig = plt.figure(figsize=(12, 7.5))
        fig.text(.06, .935, "FGO runtime breakdown", fontsize=22,
                 fontweight="bold", color="#16283E")
        fig.text(.06, .885,
                 f"Window W = {window}  |  100 epochs  |  Mean of {summary['repeats']} timed runs",
                 fontsize=12, color="#526173")
        for policy, left in (("schur", .06), ("discard", .55)):
            row = lookup[(window, policy)]
            ax = fig.add_axes([left + .035, .335, .35, .44])
            durations, percentages = pie(ax, row)
            fig.text(left + .21, .817, METHODS[policy], ha="center",
                     fontsize=17, fontweight="bold", color="#243244")
            fig.text(left + .21, .777, f"Three-stage time: {durations.sum():,.2f} ms",
                     ha="center", fontsize=11, color="#526173")
            labels = ("Add state / factor", "Estimate",
                      "Marginalize" if policy == "schur" else "Direct discard")
            for index, label in enumerate(labels):
                y = .288 - index * .049
                fig.text(left, y, "■", color=COLORS[index], fontsize=14)
                fig.text(left + .027, y, label, color="#243244", fontsize=11)
                fig.text(left + .32, y, f"{durations[index]:,.2f} ms", ha="right",
                         color="#243244", fontsize=11)
                fig.text(left + .42, y, percent(percentages[index], 2), ha="right",
                         color="#243244", fontsize=11, fontweight="bold")
        overhead = [lookup[(window, policy)]["other_mean_ms"] for policy in METHODS]
        fig.text(.06, .108,
                 f"Percentages use only the three measured stages. Unassigned overhead: "
                 f"Schur {overhead[0]:.2f} ms; discard {overhead[1]:.2f} ms.",
                 fontsize=9.1, color="#526173")
        fig.text(.06, .069,
                 "Add includes prediction and factor construction. Internal Schur work is counted only under Marginalize.",
                 fontsize=9.1, color="#526173")
        fig.text(.06, .030,
                 settings,
                 fontsize=9.1, color="#526173")
        save(fig, f"stage_timing_w{window}")

    for policy in METHODS:
        columns = min(4, len(windows))
        rows = (len(windows) + columns - 1) // columns
        fig, axes = plt.subplots(rows, columns, figsize=(3.8 * columns, 3.4 * rows + 1.7),
                                 squeeze=False)
        fig.subplots_adjust(left=.04, right=.97, bottom=.12, top=.80, hspace=.40, wspace=.12)
        fig.text(.05, .94, f"Runtime breakdown — {METHODS[policy]}", fontsize=21,
                 fontweight="bold", color="#16283E")
        fig.text(.05, .89, f"100 epochs per window · Mean of {summary['repeats']} runs · Three-stage time shares",
                 fontsize=12, color="#526173")
        labels = ["Add state / factor", "Estimate", "Marginalize" if policy == "schur" else "Direct discard"]
        fig.legend(handles=[Patch(facecolor=color, label=label) for color, label in zip(COLORS, labels)],
                   loc="upper left", bbox_to_anchor=(.04, .865), ncol=3, frameon=False)
        for ax, window in zip(axes.flat, windows):
            row = lookup[(window, policy)]
            durations, percentages = pie(ax, row, threshold=7, fontsize=12)
            ax.set_title(f"W = {window}   |   {durations.sum():,.0f} ms", fontsize=12,
                         color="#243244", pad=7)
            ax.text(.5, -.05,
                    f"Add {percent(percentages[0])}  ·  Est {percent(percentages[1])}  ·  "
                    f"{'Marg' if policy == 'schur' else 'Drop'} {percent(percentages[2])}",
                    transform=ax.transAxes, ha="center", fontsize=9, color="#526173")
        for ax in list(axes.flat)[len(windows):]:
            ax.set_visible(False)
        fig.text(.05, .055, settings + ". Unassigned overhead excluded.",
                 fontsize=10, color="#526173")
        fig.text(.05, .025, behavior,
                 fontsize=10, color="#526173")
        save(fig, f"stage_timing_overview_{policy}")
