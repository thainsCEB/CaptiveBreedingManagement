#!/usr/bin/env python3

import argparse, numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Circle

def load_table(path: str) -> pd.DataFrame:
    try:
        if path.lower().endswith(".csv"):
            return pd.read_csv(path)
        else:
            return pd.read_csv(path, sep=r"\s+|\t|,", engine="python")
    except Exception as e:
        raise SystemExit(f"Failed to read table '{path}': {e}")

def make_layout(width, height, legend_pos="right", hjust=0.22, group_bars="none"):
    rows = 1 if group_bars == "none" else 2
    if legend_pos == "none":
        gs = GridSpec(rows, 1, height_ratios=[1,0.12] if rows==2 else [1], hspace=0.02)
        fig = plt.figure(figsize=(width, height))
        ax  = fig.add_subplot(gs[0,0])
        gax = fig.add_subplot(gs[1,0]) if rows==2 else None
        return fig, ax, None, gax

    if legend_pos in ("left","right"):
        try:
            hjust = float(hjust)
        except Exception:
            hjust = 0.22
        hjust = max(0.05, min(hjust, 0.45))
        gs = GridSpec(rows, 2,
                      width_ratios=[1-hjust, hjust] if legend_pos=="right" else [hjust, 1-hjust],
                      height_ratios=[1, 0.12] if rows==2 else [1],
                      wspace=0.03, hspace=0.02)
        fig = plt.figure(figsize=(width, height))
        if legend_pos=="left":
            lax = fig.add_subplot(gs[0,0]); lax.axis("off")
            ax  = fig.add_subplot(gs[0,1])
            gax = fig.add_subplot(gs[1,1]) if rows==2 else None
        else:
            ax  = fig.add_subplot(gs[0,0])
            lax = fig.add_subplot(gs[0,1]); lax.axis("off")
            gax = fig.add_subplot(gs[1,0]) if rows==2 else None
        return fig, ax, lax, gax

    # top/bottom
    base = [0.18, 0.82] if legend_pos=="top" else [0.82, 0.18]
    if rows==2: base = base + [0.12]
    gs = GridSpec(len(base), 1, height_ratios=base, hspace=0.03)
    fig = plt.figure(figsize=(width, height))
    if legend_pos=="top":
        lax = fig.add_subplot(gs[0,0]); lax.axis("off")
        ax  = fig.add_subplot(gs[1,0])
        gax = fig.add_subplot(gs[2,0]) if rows==2 else None
    else:
        ax  = fig.add_subplot(gs[0,0])
        lax = fig.add_subplot(gs[1,0]); lax.axis("off")
        gax = fig.add_subplot(gs[2,0]) if rows==2 else None
    return fig, ax, lax, gax

def place_legend(ax, lax, legend_pos, ncols=1):
    if legend_pos == "none": return None
    handles, labels = ax.get_legend_handles_labels()
    return lax.legend(handles, labels, ncol=ncols, loc="center")

def split_by_source(df, ref_label="Reference", target_label="Target", ref_side="left", order_mode="clean", q_cols=None):
    left  = df[df["Source"].astype(str).str.lower() == ref_label.lower()].copy()
    right = df[df["Source"].astype(str).str.lower() == target_label.lower()].copy()
    other = df[~df.index.isin(left.index) & ~df.index.isin(right.index)].copy()

    def order_block(block):
        if order_mode != "clean" or q_cols is None or len(block)==0:
            return block
        q = block[q_cols].to_numpy()
        dom_idx = q.argmax(axis=1)
        dom_val = q.max(axis=1)
        block = block.copy()
        block["_dom_idx"] = dom_idx
        block["_dom_val"] = dom_val
        block.sort_values(by=["Group","_dom_idx","_dom_val","ID"], ascending=[True, True, False, True], inplace=True)
        block.drop(columns=["_dom_idx","_dom_val"], inplace=True)
        return block

    left  = order_block(left)
    right = order_block(right)
    other = order_block(other)

    if ref_side == "left":
        ordered = pd.concat([left, right, other], axis=0).reset_index(drop=True)
        ref_start = 0; ref_end = len(left)
    else:
        ordered = pd.concat([right, left, other], axis=0).reset_index(drop=True)
        ref_start = 0; ref_end = len(right)
    ordered["x"] = range(len(ordered))
    return ordered, ref_start, ref_end

def draw_connector_labels(ax, ordered, by_col, include_mask,
                          tiers=3, base_y=-0.06, delta_y=0.035,
                          pad=0.45, label_offset=0.015, linewidth=1.2, fontsize=None):
    """One connector per contiguous run of value(by_col) across included rows, staggered to avoid overlap."""
    import matplotlib.transforms as mtransforms
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

    # Build runs left to right using include_mask
    seq = [(i, ordered.iloc[i][by_col] if include_mask[i] else None) for i in range(len(ordered))]
    runs = []
    i = 0
    while i < len(seq):
        x, v = seq[i]
        if v is None:
            i += 1; continue
        x0 = x; value = v; j = i+1
        while j < len(seq) and seq[j][1] == value:
            j += 1
        x1 = seq[j-1][0]
        runs.append((x0, x1, value))
        i = j

    # Stagger tiers round-robin
    min_y = 0.0
    for r, (x0, x1, val) in enumerate(runs):
        y = base_y - (r % max(1, tiers)) * delta_y
        min_y = min(min_y, y - label_offset - 0.02)
        ax.plot([x0 - pad, x1 + pad], [y, y], transform=trans, linewidth=linewidth)
        ax.text((x0+x1)/2.0, y - label_offset, str(val), transform=trans, ha="center", va="top", fontsize=fontsize)
    return min_y

def main():
    ap = argparse.ArgumentParser(description="Ancestry (Q-matrix) barplot with clean ordering & connector labels.")
    ap.add_argument("-t","--table", required=True, help="Input (.tsv/.csv). Columns: ID, Group, Source[, Subsource], then Q columns.")
    ap.add_argument("-o","--out", required=True, help="Output image (.png/.pdf/.svg).")

    # Layout & style
    ap.add_argument("--width", type=float, default=16.0)
    ap.add_argument("--height", type=float, default=5.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--title", default="none")
    ap.add_argument("--title-pad", type=float, default=10.0)
    ap.add_argument("--x-title", default="none")
    ap.add_argument("--y-title", default="Ancestry proportion")
    ap.add_argument("--font-size", type=float, default=11.0)
    ap.add_argument("--font-family", default=None)

    ap.add_argument("--legend", choices=["left","right","top","bottom","none"], default="right")
    ap.add_argument("--hjust", type=float, default=0.22)
    ap.add_argument("--group-bars", choices=["none","top","bottom"], default="none")
    ap.add_argument("--group-bar-labels", action="store_true")

    # Source split
    ap.add_argument("--ref-side", choices=["left","right"], default="left")
    ap.add_argument("--ref-label", default="Reference")
    ap.add_argument("--target-label", default="Target")

    # Ordering
    ap.add_argument("--order", choices=["clean","none"], default="clean")

    # Labels
    ap.add_argument("--label", choices=["none","id","group","source","subsource","both"], default="id")
    ap.add_argument("--label-rotate", type=float, default=90.0)
    ap.add_argument("--label-size", type=float, default=8.0)
    ap.add_argument("--tick-pad", type=float, default=2.0)
    ap.add_argument("--x-tick-marks", choices=["on","off"], default="off")
    ap.add_argument("--label-scope", choices=["all","highlighted","none"], default="all")
    ap.add_argument("--label-set", choices=["both","reference","target"], default="both")

    # Connector styling / staggering
    ap.add_argument("--connector-linewidth", type=float, default=1.2)
    ap.add_argument("--connector-base-mode", choices=["auto","fixed"], default="auto")
    ap.add_argument("--target-display", default="", help="If set and labeling by Source, replace the display of the Target source with this string.")
    ap.add_argument("--connector-tiers", type=int, default=3)
    ap.add_argument("--connector-base-y", type=float, default=-0.06)
    ap.add_argument("--connector-delta-y", type=float, default=0.035)
    ap.add_argument("--connector-pad", type=float, default=0.45)
    ap.add_argument("--connector-label-offset", type=float, default=0.015)

    # Highlighting
    ap.add_argument("--highlight-samples", default=None)
    ap.add_argument("--highlight-groups", default=None)
    ap.add_argument("--highlight-mode", choices=["block","outline","circle"], default="block")
    ap.add_argument("--highlight-alpha", type=float, default=0.15)
    ap.add_argument("--highlight-edgewidth", type=float, default=1.5)

    # Separators
    ap.add_argument("--subspecies-linewidth", type=float, default=1.0)
    ap.add_argument("--block-sep-linewidth", type=float, default=1.5)

    args = ap.parse_args()

    # Fonts
    mpl.rcParams["font.size"] = args.font_size
    if args.font_family:
        mpl.rcParams["font.family"] = args.font_family

    # Data
    df = load_table(args.table)
    for c in ("ID","Group","Source"):
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")
    has_subsource = "Subsource" in df.columns

    meta_cols = ["ID","Group","Source"] + (["Subsource"] if has_subsource else [])
    q_cols = [c for c in df.columns if c not in meta_cols]
    if len(q_cols) < 2:
        raise SystemExit("Need at least two Q columns beyond metadata.")
    for c in q_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    # Order and split
    ordered, ref_start, ref_end = split_by_source(df, args.ref_label, args.target_label, args.ref_side, args.order, q_cols)

    # Layout
    fig, ax, lax, gax = make_layout(args.width, args.height, legend_pos=args.legend, hjust=args.hjust, group_bars=args.group_bars)

    # Bars
    bottom = np.zeros(len(ordered))
    for c in q_cols:
        ax.bar(ordered["x"], ordered[c], bottom=bottom, label=c)
        bottom += ordered[c].to_numpy()

    # Separators
    if ref_end > 0 and ref_end < len(ordered):
        ax.axvline(ref_end - 0.5, color="k", linewidth=args.block_sep_linewidth)
    if args.subspecies_linewidth > 0 and ref_end > 1:
        prev = ordered.iloc[0]["Group"] if ref_end>0 else None
        for i in range(1, ref_end):
            g = ordered.iloc[i]["Group"]
            if g != prev:
                ax.axvline(i - 0.5, color="k", linewidth=args.subspecies_linewidth, linestyle="--", alpha=0.8)
            prev = g

    # Titles
    ax.set_ylabel("" if args.y_title.lower()=="none" else args.y_title)
    if args.x_title.lower() != "none":
        ax.set_xlabel(args.x_title)
    if args.title and args.title.lower() != "none":
        ax.set_title(args.title, pad=args.title_pad)

    # Highlights
    def parse_list(arg):
        if arg is None or str(arg).strip() == "": return set()
        return set(x.strip() for x in str(arg).split(",") if x.strip())
    hi_ids = parse_list(args.highlight_samples)
    hi_grps = parse_list(args.highlight_groups)
    idxs = [i for i, r in ordered.iterrows() if (r["ID"] in hi_ids) or (r["Group"] in hi_grps)]

    
# Labeling
    if args.label == "none":
        ax.set_xticks([])
        ax.tick_params(axis="x", which="both", length=0, bottom=False, top=False)

    elif args.label in ("group","source","subsource"):
        # connector labels with staggering + label-set filtering
        ax.set_xticks([])
        ax.tick_params(axis="x", which="both", length=0, bottom=False, top=False)

        def source_ok(s):
            s = str(s).lower()
            if args.label_set == "both": return True
            if args.label_set == "reference": return s == args.ref_label.lower()
            if args.label_set == "target": return s == args.target_label.lower()
            return True

        include_mask = []
        for i, row in ordered.iterrows():
            in_scope = True if args.label_scope != "highlighted" else (i in idxs)
            if args.label_scope == "none":
                in_scope = False
            include_mask.append(in_scope and source_ok(row["Source"]))

        display_col = "Group" if args.label=="group" else ("Subsource" if args.label=="subsource" else "Source")
        ordered_disp = ordered.copy()
        if display_col == "Source" and isinstance(args.target_display, str) and args.target_display.strip():
            mask_target = ordered_disp["Source"].astype(str).str.lower() == args.target_label.lower()
            ordered_disp.loc[mask_target, "Source"] = args.target_display

        # Auto baseline hugging the axis unless fixed is requested
        ymin, ymax = ax.get_ylim()
        if args.connector_base_mode == "auto":
            auto_base = -0.012 * max(8.0, float(args.label_size)) / 8.0
            base_y = auto_base
        else:
            base_y = args.connector_base_y

        if ymin >= 0:
            ax.set_ylim(min(base_y - (args.connector_tiers * args.connector_delta_y) - 0.06, ymin), max(1.0001, ymax))

        min_y = draw_connector_labels(ax, ordered_disp, by_col=display_col, include_mask=include_mask,
                                      tiers=args.connector_tiers, base_y=base_y, delta_y=args.connector_delta_y,
                                      pad=args.connector_pad, label_offset=args.connector_label_offset,
                                      linewidth=args.connector_linewidth, fontsize=max(8, int(args.label_size)))
        ymin2, ymax2 = ax.get_ylim()
        if ymin2 > min_y:
            ax.set_ylim(min_y, ymax2)

    else:
        # per-bar labels with scope & set
        def source_ok(s):
            s = str(s).lower()
            if args.label_set == "both": return True
            if args.label_set == "reference": return s == args.ref_label.lower()
            if args.label_set == "target": return s == args.target_label.lower()
            return True
        scope = args.label_scope
        if scope not in ("all","highlighted","none"): scope = "all"
        mask = [True]*len(ordered) if scope=="all" else ([i in idxs for i in ordered.index] if scope=="highlighted" else [False]*len(ordered))
        labs = []
        for i, row in ordered.iterrows():
            if not mask[i] or not source_ok(row["Source"]):
                labs.append(""); continue
            if args.label == "id":
                labs.append(str(row["ID"]))
            elif args.label == "both":
                labs.append(f"{row['Group']}:{row['ID']}")
            else:
                labs.append(str(row["ID"]))
        ax.set_xticks(range(len(ordered)))
        ax.set_xticklabels(labs, rotation=args.label_rotate, ha="right", fontsize=args.label_size)
        ax.tick_params(axis="x", pad=args.tick_pad)
        if args.x_tick_marks == "off":
            ax.tick_params(axis="x", which="both", length=0, bottom=False, top=False)

# Optional group bars panel
    if args.group_bars in ("top","bottom") and gax is not None:
        gax.axis("off")
        for grp, df_g in ordered.groupby("Group"):
            xs = list(df_g["x"].sort_values())
            if not xs: continue
            start = xs[0]; prev = xs[0]; runs = []
            for x in xs[1:]:
                if x == prev + 1: prev = x
                else: runs.append((start, prev)); start = x; prev = x
            runs.append((start, prev))
            for (x0,x1) in runs:
                gax.add_patch(Rectangle((x0-0.5,0), (x1-x0+1), 1, linewidth=0, facecolor="0.9"))
                if args.group_bar_labels:
                    gax.text((x0+x1)/2.0, 0.5, str(grp), ha="center", va="center", fontsize=max(8,int(args.label_size)))
        gax.set_xlim(-0.5, max(ordered["x"])+0.5); gax.set_ylim(0,1)

    # Y limits & highlights
    ax.set_ylim(0, 1.0001)
    if idxs:
        if args.highlight_mode == "block":
            for i in idxs: ax.add_patch(Rectangle((i-0.5,0), 1.0, 1.0, fill=True, alpha=args.highlight_alpha, edgecolor="none"))
        elif args.highlight_mode == "outline":
            for i in idxs: ax.add_patch(Rectangle((i-0.5,0), 1.0, 1.0, fill=False, linewidth=args.highlight_edgewidth))
        elif args.highlight_mode == "circle":
            for i in idxs: ax.add_patch(Circle((i,0.5), 0.35, fill=False, linewidth=args.highlight_edgewidth))

    # Legend
    extra = None
    if args.legend != "none":
        extra = place_legend(ax, lax, args.legend, ncols=max(1, int(len(q_cols)//5)+1))

    # No outer box
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(left=False, bottom=False)

    # Save safely
    save_kwargs = dict(dpi=args.dpi, bbox_inches="tight", pad_inches=0.5, facecolor="white")
    if extra is not None:
        save_kwargs["bbox_extra_artists"] = (extra,)
    ext = args.out.lower().split(".")[-1]
    if ext in ("png","jpg","jpeg","tif","tiff","pdf","svg"):
        plt.savefig(args.out, **save_kwargs)
    else:
        plt.savefig(args.out + ".png", **save_kwargs)

if __name__ == "__main__":
    main()
