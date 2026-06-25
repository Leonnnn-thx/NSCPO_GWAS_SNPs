#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import itertools
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.ticker import EngFormatter
from matplotlib.patches import Rectangle
import cooler
import bioframe

plt.rcParams["pdf.compression"] = 9
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 11


def parse_region(region_str):
    chrom, pos = region_str.split(":")
    s, e = pos.split("-")
    start = int(s.replace(",", ""))
    end = int(e.replace(",", ""))
    return chrom, start, end


def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    n = matrix_c.shape[0]
    start_pos_vector = [start + resolution * i for i in range(n + 1)]
    t = np.array([[1, 0.5], [-1, 0.5]])

    matrix_a = np.dot(
        np.array([(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)]), t
    )
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)

    im = ax.pcolormesh(x, y, np.flipud(matrix_c), shading="auto", *args, **kwargs)
    im.set_rasterized(True)
    return im


bp_formatter = EngFormatter("b")


def format_ticks(ax, x=True, y=True, rotate=False):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis="x", rotation=45)


def extract_tads(ins_table, is_boundary_col, max_tad_length=3_000_000):
    need = ["chrom", "start", "end", is_boundary_col]
    for c in need:
        if c not in ins_table.columns:
            raise ValueError(f"缺少列: {c}")

    tads = bioframe.merge(ins_table.loc[ins_table[is_boundary_col] == False, ["chrom", "start", "end"]])
    tads = tads[(tads["end"] - tads["start"]) <= max_tad_length].reset_index(drop=True)
    return tads[["chrom", "start", "end"]]


def build_tad_overlay_matrix(region_start, region_end, resolution, tad_df):
    n_bins = int(np.ceil((region_end - region_start) / resolution))
    mat = np.full((n_bins, n_bins), np.nan, dtype=float)

    for _, row in tad_df.iterrows():
        s = max(int(row["start"]), region_start)
        e = min(int(row["end"]), region_end)
        if e <= s:
            continue

        i0 = int((s - region_start) // resolution)
        i1 = int(np.ceil((e - region_start) / resolution))
        i0 = max(i0, 0)
        i1 = min(i1, n_bins)

        if i1 - i0 < 2:
            continue

        mat[i0:i1, i0:i1] = 1.0
        if i1 - i0 > 2:
            mat[i0 + 1:i1 - 1, i0 + 1:i1 - 1] = np.nan

    return mat


def parse_gtf_attributes(attr_str):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)";', attr_str):
        d[m.group(1)] = m.group(2)
    return d


def load_gtf_features_in_region(gtf_path, chrom, start, end):
    """
    返回:
      genes_df: gene/transcript区间（优先gene；若缺gene可退化用transcript）
      exons_df: exon区间
    """
    colnames = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    usecols = [0, 2, 3, 4, 6, 8]

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=colnames,
        usecols=usecols,
        dtype={0: str, 2: str, 3: int, 4: int, 6: str, 8: str}
    )

    cands = {chrom}
    if chrom.startswith("chr"):
        cands.add(chrom.replace("chr", ""))
    else:
        cands.add("chr" + chrom)

    gtf = gtf[gtf["seqname"].isin(cands)].copy()
    gtf = gtf[(gtf["end"] >= start) & (gtf["start"] <= end)].copy()
    if gtf.empty:
        return pd.DataFrame(), pd.DataFrame()

    attrs = gtf["attribute"].apply(parse_gtf_attributes)
    gtf["gene_id"] = attrs.apply(lambda x: x.get("gene_id", "NA"))
    gtf["gene_name"] = attrs.apply(lambda x: x.get("gene_name", x.get("gene_id", "NA")))
    gtf["transcript_id"] = attrs.apply(lambda x: x.get("transcript_id", "NA"))

    genes = gtf[gtf["feature"] == "gene"].copy()
    transcripts = gtf[gtf["feature"] == "transcript"].copy()
    exons = gtf[gtf["feature"] == "exon"].copy()

    if genes.empty and (not transcripts.empty):
        genes = transcripts.copy()
        genes["feature"] = "gene"

    if not genes.empty:
        genes = genes[["seqname", "start", "end", "strand", "gene_id", "gene_name"]].drop_duplicates()

    if not exons.empty:
        exons = exons[["seqname", "start", "end", "strand", "gene_id", "gene_name", "transcript_id"]].drop_duplicates()

    return genes, exons


def assign_gene_rows(genes_df, row_gap_bp=0):
    if genes_df.empty:
        return genes_df

    genes = genes_df.sort_values(["start", "end"]).copy()
    row_ends = []
    rows = []

    for _, r in genes.iterrows():
        placed = False
        for i in range(len(row_ends)):
            if r["start"] > row_ends[i] + row_gap_bp:
                rows.append(i)
                row_ends[i] = r["end"]
                placed = True
                break
        if not placed:
            rows.append(len(row_ends))
            row_ends.append(r["end"])

    genes["row"] = rows
    return genes


def _stagger_text_y(placed_x, x, base_y, dy=0.22, max_levels=4, min_dx=120000):
    level = 0
    while level < max_levels:
        conflict = any(abs(x - px) < min_dx and pl == level for px, pl in placed_x)
        if not conflict:
            placed_x.append((x, level))
            return base_y + level * dy
        level += 1
    placed_x.append((x, max_levels - 1))
    return base_y + (max_levels - 1) * dy


def draw_gene_track(ax, genes_df, exons_df, region_start, region_end, max_rows=20):
    ax.set_xlim(region_start, region_end)

    if genes_df.empty:
        ax.text(0.5, 0.5, "No genes in region", transform=ax.transAxes, ha="center", va="center")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return

    g = genes_df.copy()
    g["len"] = g["end"] - g["start"]
    g = g.sort_values(["gene_id", "len"], ascending=[True, False]).drop_duplicates("gene_id", keep="first")
    g = g.drop(columns=["len"]).reset_index(drop=True)

    genes_plot = assign_gene_rows(g, row_gap_bp=0)
    genes_plot = genes_plot[genes_plot["row"] < max_rows].copy()

    ex = exons_df.copy()
    if not ex.empty:
        ex = ex.drop_duplicates(subset=["gene_id", "start", "end"])

    h_exon = 0.34
    label_cache_by_row = {}

    for _, gene in genes_plot.iterrows():
        y = float(gene["row"])
        gs = max(int(gene["start"]), region_start)
        ge = min(int(gene["end"]), region_end)
        if ge <= gs:
            continue

        ax.hlines(y, gs, ge, color="#444444", linewidth=1.2, zorder=1)

        n_arrows = max(1, int((ge - gs) // 140000))
        xs = np.linspace(gs, ge, n_arrows + 2)[1:-1]
        for xx in xs:
            if gene["strand"] == "+":
                ax.plot([xx - 9000, xx, xx - 9000], [y - 0.06, y, y + 0.06], color="#444444", lw=0.8)
            elif gene["strand"] == "-":
                ax.plot([xx + 9000, xx, xx + 9000], [y - 0.06, y, y + 0.06], color="#444444", lw=0.8)

        ex_sub = ex[ex["gene_id"] == gene["gene_id"]] if not ex.empty else pd.DataFrame()
        for _, e in ex_sub.iterrows():
            es = max(int(e["start"]), region_start)
            ee = min(int(e["end"]), region_end)
            if ee <= es:
                continue
            ax.add_patch(
                Rectangle((es, y - h_exon / 2), ee - es, h_exon,
                          facecolor="#1f77b4", edgecolor="none", alpha=0.9, zorder=2)
            )

        if y not in label_cache_by_row:
            label_cache_by_row[y] = []
        y_text = _stagger_text_y(
            label_cache_by_row[y], gs, y + 0.30,
            dy=0.20, max_levels=4, min_dx=120000
        )
        ax.text(gs, y_text, str(gene["gene_name"]), fontsize=8, ha="left", va="bottom", color="#222222")

    ymax = max(1, int(genes_plot["row"].max()) + 1 if not genes_plot.empty else 1)
    ax.set_ylim(-0.8, ymax + 1.2)
    ax.set_ylabel("Genes")
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def parse_kv_map(items, value_split=",", cast_key=int, cast_val=int):
    """
    e.g.
      value_split=None:
        ["25600:/path/a.tsv", "12800:/path/b.tsv"] -> {25600: "/path/a.tsv", ...}
      value_split=",":
        ["25600:256000,640000"] -> {25600:[256000,640000]}
    """
    out = {}
    for it in items:
        if ":" not in it:
            raise ValueError(f"格式错误（缺少:）: {it}")
        k, v = it.split(":", 1)
        k = cast_key(k)
        if value_split is None:
            out[k] = v
        else:
            vals = [x for x in v.split(value_split) if x.strip() != ""]
            out[k] = [cast_val(x) for x in vals]
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Hi-C + multi-resolution TAD score tracks + TAD overlay + gene track"
    )
    parser.add_argument("-i", "--input", required=True, help="Input mcool file")
    parser.add_argument("-r", "--region", required=True, help="chr:start-end")
    parser.add_argument("--gtf", required=True, help="Reference GTF file")
    parser.add_argument("-o", "--output-prefix", required=True, help="output prefix")

    # 新增：热图分辨率与多分辨率score
    parser.add_argument("--hic-res", type=int, required=True, help="Resolution for Hi-C heatmap, e.g. 3200")
    parser.add_argument("--score-res-list", type=int, nargs="+", required=True,
                        help="Resolutions for TAD score panels, e.g. 25600 12800 6400 3200")
    parser.add_argument("--insulation-map", nargs="+", required=True,
                        help='Mapping res:path, e.g. 25600:/a.tsv 12800:/b.tsv')
    parser.add_argument("--windows-map", nargs="+", required=True,
                        help='Mapping res:w1,w2,... e.g. 25600:256000,640000,1280000,2560000')

    parser.add_argument("--max-tad-length", type=int, default=3000000)
    parser.add_argument("--overlay-n-tads", type=int, default=20)
    parser.add_argument("--overlay-boundary-res", type=int, default=None,
                        help="Use which score resolution table to define TAD overlay boundaries; default=max(score-res)")
    parser.add_argument("--gene-max-rows", type=int, default=20)

    parser.add_argument("--vmin", type=float, default=1e-3)
    parser.add_argument("--vmax", type=float, default=1e-1)
    parser.add_argument("--cmap", type=str, default="YlOrRd")
    parser.add_argument("--dpi", type=int, default=220)

    args = parser.parse_args()

    chrom, start, end = parse_region(args.region)
    region = (chrom, start, end)
    L = end - start

    ins_map = parse_kv_map(args.insulation_map, value_split=None, cast_key=int, cast_val=str)
    win_map = parse_kv_map(args.windows_map, value_split=",", cast_key=int, cast_val=int)

    score_res_list = sorted(args.score_res_list, reverse=True)  # 高->低
    for rr in score_res_list:
        if rr not in ins_map:
            raise ValueError(f"score res={rr} 未在 --insulation-map 提供")
        if rr not in win_map:
            raise ValueError(f"score res={rr} 未在 --windows-map 提供")

    # overlay来源分辨率：默认最大res（最平滑）
    overlay_res = args.overlay_boundary_res if args.overlay_boundary_res is not None else max(score_res_list)
    if overlay_res not in ins_map:
        raise ValueError(f"overlay boundary res={overlay_res} 未在 --insulation-map 提供")
    if overlay_res not in win_map or len(win_map[overlay_res]) == 0:
        raise ValueError(f"overlay boundary res={overlay_res} 未在 --windows-map 提供有效窗口")

    # 读取Hi-C（热图）
    hic_res = args.hic_res
    uri = f"{args.input}::resolutions/{hic_res}"
    clr = cooler.Cooler(uri)
    try:
        data = clr.matrix(balance=True).fetch(region)
    except Exception:
        print("[WARN] balance=True失败，改用balance=False")
        data = clr.matrix(balance=False).fetch(region)
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    # 用 overlay_res 的 insulation table 提取TAD并覆盖在热图上
    ins_overlay = pd.read_csv(ins_map[overlay_res], sep="\t")
    ins_overlay = ins_overlay[ins_overlay["chrom"] == chrom].copy().reset_index(drop=True)
    if ins_overlay.empty:
        raise ValueError(f"overlay表在染色体 {chrom} 为空: res={overlay_res}")

    w0 = win_map[overlay_res][0]
    is_boundary_col = f"is_boundary_{w0}"
    if is_boundary_col not in ins_overlay.columns:
        raise ValueError(f"overlay表缺少列: {is_boundary_col} (res={overlay_res})")

    tads_table = extract_tads(ins_overlay, is_boundary_col, max_tad_length=args.max_tad_length)
    tads_region = bioframe.select(tads_table, region).sort_values(["start", "end"]).reset_index(drop=True)
    tads_region_n = tads_region.head(args.overlay_n_tads)
    overlay_mat = build_tad_overlay_matrix(start, end, hic_res, tads_region_n)

    # gene track
    genes_df, exons_df = load_gtf_features_in_region(args.gtf, chrom, start, end)

    # 动态布局：top热图 + N个score + gene
    n_score = len(score_res_list)
    height_ratios = [6.2] + [1.0] * n_score + [1.6]
    fig, axes = plt.subplots(
        2 + n_score, 1,
        figsize=(14, 5.0 + 1.15 * n_score),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.08}
    )

    ax_top = axes[0]
    score_axes = axes[1:1 + n_score]
    ax_gene = axes[-1]

    norm = LogNorm(vmin=args.vmin, vmax=args.vmax)

    # top: Hi-C + overlay
    im = pcolormesh_45deg(ax_top, data, start=start, resolution=hic_res, norm=norm, cmap=args.cmap)
    overlay_cmap = ListedColormap(["#b7e4f9"])
    pcolormesh_45deg(
        ax_top, overlay_mat, start=start, resolution=hic_res,
        cmap=overlay_cmap, vmin=0, vmax=1, alpha=0.75
    )

    ax_top.set_xlim(start, end)
    ax_top.set_ylim(0, L)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.set_yticks([0, L / 2.0, L])
    ax_top.set_yticklabels(["0 b", f"{int(L/2/1000)} kb", f"{int(L/1000)} kb"])
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.set_title(
        f"Hi-C + TAD overlay | {chrom}:{start}-{end} | hic_res={hic_res} bp | overlay_res={overlay_res} bp",
        fontsize=12
    )

    # 多分辨率 score
    colors = ["#1f4aff", "#8a2be2", "#ff4f81", "#ff8c33", "#2ca02c", "#d62728"]

    for ax_s, rr in zip(score_axes, score_res_list):
        ins_table = pd.read_csv(ins_map[rr], sep="\t")
        ins_table = ins_table[ins_table["chrom"] == chrom].copy().reset_index(drop=True)
        insul_region = bioframe.select(ins_table, region)
        if insul_region.empty:
            raise ValueError(f"区域 {region} 在 insulation table 中为空: res={rr}")

        xmid = insul_region[["start", "end"]].mean(axis=1)
        plotted = 0

        for i, w in enumerate(win_map[rr]):
            col = f"log2_insulation_score_{w}"
            if col in insul_region.columns:
                ax_s.plot(xmid, insul_region[col], lw=1.35, color=colors[i % len(colors)], label=f"{w}")
                plotted += 1

        if plotted == 0:
            raise ValueError(f"res={rr} 未找到任何 log2_insulation_score_<window> 列，请检查 --windows-map")

        ax_s.spines["top"].set_visible(False)
        ax_s.spines["right"].set_visible(False)
        ax_s.set_ylabel(f"TAD\n{rr}bp", rotation=0, labelpad=25, va="center")
        ax_s.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_s.legend(
            loc="upper left",
            ncol=min(4, plotted),
            frameon=True,
            fontsize=8,
            title=f"res={rr}",
            title_fontsize=8
        )

    # gene
    draw_gene_track(ax_gene, genes_df, exons_df, start, end, max_rows=args.gene_max_rows)
    ax_gene.set_xlabel(f"{chrom} position (bp)")
    format_ticks(ax_gene, x=True, y=False, rotate=False)

    # layout
    fig.subplots_adjust(left=0.07, right=0.92, top=0.95, bottom=0.08, hspace=0.08)

    # colorbar
    pos = ax_top.get_position()
    cax = fig.add_axes([
        pos.x1 + 0.010,
        pos.y0 + 0.08 * pos.height,
        0.012,
        0.84 * pos.height
    ])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Contact Frequency")

    out_png = f"{args.output_prefix}.hic_multiResTAD_gene.png"
    out_pdf = f"{args.output_prefix}.hic_multiResTAD_gene.pdf"
    out_tad = f"{args.output_prefix}.region_tads.overlayRes{overlay_res}.tsv"

    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    tads_region.to_csv(out_tad, sep="\t", index=False)

    print("[DONE] 输出文件：")
    print(" -", out_png)
    print(" -", out_pdf)
    print(" -", out_tad)


if __name__ == "__main__":
    main()
