#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, Circle
from matplotlib.gridspec import GridSpec
import argparse

def load_table(path):
    try:
        if path.lower().endswith(".csv"):
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    return df

def compute_order(df, ref_label="Reference", target_label="Target", ref_side="left"):
    left = df[df["Source"].astype(str).str.lower() == ref_label.lower()].copy()
    right = df[df["Source"].astype(str).str.lower() == target_label.lower()].copy()
    left = left.sort_values(["Group","ID"])
    right = right.sort_values(["Group","ID"])
    if ref_side == "left":
        ordered = pd.concat([left, right], axis=0).reset_index(drop=True)
        ref_start = 0
        ref_end = len(left)
    else:
        ordered = pd.concat([right, left], axis=0).reset_index(drop=True)
        ref_start = len(right)
        ref_end = len(right) + len(left)
    ordered["x"] = range(len(ordered))
    return ordered, ref_start, ref_end, left, right

def group_breaks_in_span(ordered, span_start, span_end):
    if span_end - span_start <= 1:
        return []
    sub = ordered.iloc[span_start:span_end].reset_index(drop=True)
    brks = []
    for i in range(len(sub)-1):
        if sub.loc[i, "Group"] != sub.loc[i+1, "Group"]:
            brks.append(span_start + i + 0.5)
    return brks

def parse_list(arg):
    if arg is None or str(arg).strip() == "":
        return set()
    return set([x.strip() for x in str(arg).split(",") if x.strip()])

def contiguous_runs(indices):
    if not indices:
        return []
    indices = sorted(indices)
    runs = []
    start = indices[0]
    prev = indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
            continue
        runs.append((start, prev))
        start = i
        prev = i
    runs.append((start, prev))
    return runs

def make_layout(width, height, legend_pos, group_bars, hjust):
    # hjust controls legend panel width fraction for left/right legends.
    try:
        hjust = float(hjust)
    except Exception:
        hjust = 0.22
    hjust = max(0.05, min(hjust, 0.45))

    if legend_pos in ("top","bottom"):
        if group_bars == "none":
            ratios = [0.18, 0.82] if legend_pos=="top" else [0.82, 0.18]
            gs = GridSpec(2, 1, height_ratios=ratios, hspace=0.03)
            fig = plt.figure(figsize=(width, height))
            if legend_pos=="top":
                lax = fig.add_subplot(gs[0,0]); lax.axis("off")
                ax  = fig.add_subplot(gs[1,0])
            else:
                ax  = fig.add_subplot(gs[0,0])
                lax = fig.add_subplot(gs[1,0]); lax.axis("off")
            return fig, ax, lax, None
        else:
            if legend_pos=="top" and group_bars=="top":
                ratios = [0.16, 0.10, 0.74]
                order  = ("legend","group","main")
            elif legend_pos=="top" and group_bars=="bottom":
                ratios = [0.16, 0.74, 0.10]
                order  = ("legend","main","group")
            elif legend_pos=="bottom" and group_bars=="top":
                ratios = [0.10, 0.74, 0.16]
                order  = ("group","main","legend")
            else:  # bottom + bottom
                ratios = [0.74, 0.10, 0.16]
                order  = ("main","group","legend")
            gs = GridSpec(3, 1, height_ratios=ratios, hspace=0.03)
            fig = plt.figure(figsize=(width, height))
            slots = {}
            for i, name in enumerate(order):
                slots[name] = fig.add_subplot(gs[i,0])
                if name in ("legend","group"):
                    slots[name].axis("off")
            return fig, slots["main"], slots["legend"], slots["group"]
    else:
        if group_bars == "none":
            gs = GridSpec(1, 2, width_ratios=[1-hjust, hjust] if legend_pos=="right" else [hjust, 1-hjust], wspace=0.03)
            fig = plt.figure(figsize=(width, height))
            if legend_pos=="left":
                lax = fig.add_subplot(gs[0,0]); lax.axis("off")
                ax  = fig.add_subplot(gs[0,1])
            else:
                ax  = fig.add_subplot(gs[0,0])
                lax = fig.add_subplot(gs[0,1]); lax.axis("off")
            return fig, ax, lax, None
        else:
            gs = GridSpec(2, 2,
                          height_ratios=[0.12, 0.88],
                          width_ratios=[1-hjust, hjust] if legend_pos=="right" else [hjust, 1-hjust],
                          hspace=0.02, wspace=0.03)
            fig = plt.figure(figsize=(width, height))
            gax = fig.add_subplot(gs[0, :]); gax.axis("off")
            if legend_pos=="left":
                lax = fig.add_subplot(gs[1,0]); lax.axis("off")
                ax  = fig.add_subplot(gs[1,1])
            else:
                ax  = fig.add_subplot(gs[1,0])
                lax = fig.add_subplot(gs[1,1]); lax.axis("off")
            return fig, ax, lax, gax

def place_legend(fig, ax, lax, ncols, pos):
    handles, labels = ax.get_legend_handles_labels()
    if pos == "none":
        return
    if lax is None:
        fig.legend(handles, labels, ncol=ncols, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    elif pos in ("top","bottom"):
        lax.legend(handles, labels, ncol=ncols, loc="center")
    else:
        lax.legend(handles, labels, ncol=1, loc="center")

def draw_group_bars(gax, ordered, show_labels=True):
    if gax is None:
        return
    runs = []
    start = 0
    for i in range(len(ordered)-1):
        if ordered.loc[i, "Group"] != ordered.loc[i+1, "Group"]:
            runs.append((start, i, ordered.loc[i, "Group"]))
            start = i+1
    runs.append((start, len(ordered)-1, ordered.loc[len(ordered)-1, "Group"]))
    gax.set_xlim(-0.5, len(ordered)-0.5)
    gax.set_ylim(0, 1)
    for a, b, grp in runs:
        rect = Rectangle((a-0.5, 0.1), (b-a+1), 0.8, fill=False, linewidth=1.0)
        gax.add_patch(rect)
        if show_labels:
            cx = (a + b) / 2.0
            gax.text(cx, 0.5, str(grp), ha="center", va="center")
    gax.axis("off")

def main():
    ap = argparse.ArgumentParser(description="Admixture barplot with reference on one side and target on the other.")
    ap.add_argument("--table", "-t", required=True, help="Input TSV/CSV with columns: ID, Group, Source, and K columns (ancestry components).")
    ap.add_argument("--out", "-o", required=True, help="Output image (.png/.pdf/.svg).")
    ap.add_argument("--ref-side", choices=["left","right"], default="left", help="Place reference block on left or right. Default: left")
    ap.add_argument("--ref-label", default="Reference", help="Label used in Source column for reference samples. Default: Reference")
    ap.add_argument("--target-label", default="Target", help="Label used in Source column for target samples. Default: Target")
    ap.add_argument("--width", type=float, default=16.0, help="Figure width (inches). Default: 16")
    ap.add_argument("--height", type=float, default=5.0, help="Figure height (inches). Default: 5")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for raster outputs. Default: 300")
    ap.add_argument("--font-size", type=float, default=11.0, help="Base font size. Default: 11")
    ap.add_argument("--font-family", default=None, help="Font family (e.g., DejaVu Sans, Arial). Optional.")
    ap.add_argument("--legend", choices=["top","bottom","left","right","best","none"], default="right", help="Legend placement. Default: right")
    ap.add_argument("--hjust", type=float, default=0.22, help="When legend is left/right, fraction of figure width reserved for the legend panel (distance from axis). Range ~0.05-0.45. Default: 0.22")
    # Titles and padding (defaults: no main title, no x-title)
    ap.add_argument("--title", default="none", help="Main title (default none).")
    ap.add_argument("--title-pad", type=float, default=8.0, help="Padding above the title. Default: 8")
    ap.add_argument("--x-title", default="none", help="X-axis title (default none).")
    ap.add_argument("--x-title-pad", type=float, default=8.0, help="X-axis title padding. Default: 8")
    ap.add_argument("--y-title", default="Ancestry proportion", help="Y-axis title. Use 'none' to hide.")
    ap.add_argument("--y-title-pad", type=float, default=6.0, help="Y-axis title padding. Default: 6")
    # Labels
    ap.add_argument("--label", choices=["none","id","group","both"], default="id", help="Per-sample label content. Default: id")
    ap.add_argument("--label-scope", choices=["all","highlighted","none"], default="all", help="Which samples to label. Default: all")
    ap.add_argument("--label-set", choices=["both","reference","target"], default="both", help="Restrict which samples (by Source) can display labels. Default: both")
    ap.add_argument("--label-rotate", type=float, default=90.0, help="Rotation angle for labels if shown. Default: 90")
    ap.add_argument("--label-size", type=float, default=8.0, help="Font size for labels if shown. Default: 8")
    ap.add_argument("--tick-pad", type=float, default=2.0, help="Padding between tick labels and axis. Default: 2")
    ap.add_argument("--x-tick-marks", choices=["on","off"], default="off", help="Show small tick marks on the x-axis. Default: off")
    # Highlighting
    ap.add_argument("--highlight-samples", default=None, help="Comma-separated sample IDs to highlight (and optionally label).")
    ap.add_argument("--highlight-groups", default=None, help="Comma-separated group names to highlight (and optionally label).")
    ap.add_argument("--highlight-mode", choices=["block","outline","circle"], default="block", help="Highlight style. Default: block")
    ap.add_argument("--highlight-alpha", type=float, default=0.15, help="Alpha for block highlights. Default: 0.15")
    ap.add_argument("--highlight-edgewidth", type=float, default=1.5, help="Edge width for outline/circle. Default: 1.5")
    # Group bars
    ap.add_argument("--group-bars", choices=["none","top","bottom"], default="none", help="Add horizontal bars to identify groups. Default: none")
    ap.add_argument("--group-bar-labels", action="store_true", help="Show text labels within group bars.")
    # Separator line widths
    ap.add_argument("--subspecies-linewidth", type=float, default=1.0, help="Line width for dashed separators within the reference block. Default: 1.0")
    ap.add_argument("--block-sep-linewidth", type=float, default=1.5, help="Line width for the solid separator between reference and target. Default: 1.5")
    args = ap.parse_args()

    # Fonts
    mpl.rcParams["font.size"] = args.font_size
    if args.font_family:
        mpl.rcParams["font.family"] = args.font_family

    df = load_table(args.table)

    required = ["ID","Group","Source"]
    for col in required:
        if col not in df.columns:
            raise SystemExit(f"Missing required column: {col}")
    anc_cols = [c for c in df.columns if c not in required]
    if len(anc_cols) < 2:
        raise SystemExit("Need at least two ancestry columns after ID, Group, Source.")
    for c in anc_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    ordered, ref_start, ref_end, left_df, right_df = compute_order(df, args.ref_label, args.target_label, args.ref_side)

    # Create layout with dedicated legend and optional group bars axes
    fig, ax, lax, gax = make_layout(args.width, args.height, args.legend, args.group_bars, args.hjust)
    ax.grid(False)

    # Stacked bars
    bottom_vals = pd.Series([0.0]*len(ordered), index=ordered.index)
    for c in anc_cols:
        ax.bar(ordered["x"], ordered[c], bottom=bottom_vals, label=c)
        bottom_vals = bottom_vals + ordered[c]

    # Dashed lines inside reference only
    dashes = group_breaks_in_span(ordered, ref_start, ref_end)
    import matplotlib.transforms as mtransforms
    _trans_v = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for b in dashes:
        ax.plot([b, b], [-0.04, 1.04], transform=_trans_v,
                linewidth=args.subspecies_linewidth, linestyle=(0,(4,4)), color="black", clip_on=False)

    # Solid separator between reference and target
    if 0 < ref_end < len(ordered):
        import matplotlib.transforms as mtransforms
        _trans_v = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        _x = ref_end - 0.5
        ax.plot([_x, _x], [-0.04, 1.04], transform=_trans_v,
                color="black", linewidth=args.block_sep_linewidth, clip_on=False)

    # Highlights
    h_samples = parse_list(args.highlight_samples)
    h_groups  = parse_list(args.highlight_groups)
    idxs = []
    if h_samples:
        idxs.extend(ordered.index[ordered["ID"].astype(str).isin(h_samples)].tolist())
    if h_groups:
        idxs.extend(ordered.index[ordered["Group"].astype(str).isin(h_groups)].tolist())
    idxs = sorted(set(idxs))

    if idxs:
        if args.highlight_mode == "block":
            for a, b in contiguous_runs(idxs):
                rect = Rectangle((a-0.5, 0), (b-a+1), 1.0, facecolor="0.3", alpha=args.highlight_alpha, edgecolor=None, linewidth=0)
                ax.add_patch(rect)
        elif args.highlight_mode == "outline":
            for i in idxs:
                rect = Rectangle((i-0.5, 0), 1.0, 1.0, fill=False, edgecolor="black", linewidth=args.highlight_edgewidth)
                ax.add_patch(rect)
        elif args.highlight_mode == "circle":
            for i in idxs:
                circ = Circle((i, 0.5), 0.48, fill=False, edgecolor="black", linewidth=args.highlight_edgewidth)
                ax.add_patch(circ)

    # Axis titles with padding and option to hide
    if args.y_title and args.y_title.lower() != "none":
        ax.set_ylabel(args.y_title, labelpad=args.y_title_pad)
    if args.x_title and args.x_title.lower() != "none":
        ax.set_xlabel(args.x_title, labelpad=args.x_title_pad)

    if args.title and args.title.lower() != "none":
        ax.set_title(args.title, pad=args.title_pad)

    # Legend
    place_legend(fig, ax, lax, ncols=len(anc_cols), pos=args.legend)

    # Optional group bars
    draw_group_bars(gax, ordered, show_labels=args.group_bar_labels)

    # X tick labels
    show_labels = args.label != "none" and args.label_scope != "none"
    if show_labels:
        # Base mask from label-scope
        if args.label_scope == "highlighted" and len(idxs) > 0:
            scope_mask = [i in idxs for i in ordered.index]
        else:
            scope_mask = [True]*len(ordered)
        # Source-based mask from label-set
        def allowed_source(src):
            s = str(src).lower()
            if args.label_set == "both":
                return True
            if args.label_set == "reference":
                return s == args.ref_label.lower()
            if args.label_set == "target":
                return s == args.target_label.lower()
            return True
        # Build labels
        labs = []
        for i, row in ordered.iterrows():
            if not scope_mask[i] or not allowed_source(row["Source"]):
                labs.append("")
                continue
            if args.label == "id":
                labs.append(str(row["ID"]))
            elif args.label == "group":
                labs.append(str(row["Group"]))
            elif args.label == "both":
                labs.append(f"{row['Group']}:{row['ID']}")
            else:
                labs.append("")
        ax.set_xticks(ordered["x"])
        ax.set_xticklabels(labs, rotation=args.label_rotate, fontsize=args.label_size, ha="right")
        ax.tick_params(axis="x", pad=args.tick_pad)
        # Hide or show tick marks (stubs)
        if args.x_tick_marks == "off":
            ax.tick_params(axis="x", which="both", length=0, bottom=False, top=False)
        else:
            ax.tick_params(axis="x", which="both", bottom=True, top=False)
    else:
        ax.set_xticks([])
        ax.tick_params(axis="x", which="both", length=0, bottom=False, top=False)

    ax.set_ylim(0, 1.0001)

    # Save
    ext = args.out.lower().split(".")[-1]
    if ext in ("png","jpg","jpeg","tif","tiff","pdf","svg"):
        plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    else:
        plt.savefig(args.out + ".png", dpi=args.dpi, bbox_inches="tight")

if __name__ == "__main__":
    main()
