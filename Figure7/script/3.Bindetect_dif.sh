#!/bin/bash
# bindetect_comparison.sh - Flexible TOBIAS BINDetect comparison script
"""
example usage:
bash Bindetect_dif.sh -1 HEPM_conEnh_footprint_30top_step4_merged.bed 
                    -2 HEPM_conEnhancers_Footprint_step4_merged_bottom70p.bed 
                    -f HEPM_ATAC_merged_footprint_con.bw -m /home/thx/NSCPO_Variant/Hierarchy_Enh_net/HEPM_expressed_TFs.jaspar 
                    -g /home/thx/reference_genome/hg19/hg19.fa 
                    -o ConEnh 
                    -p conEnh_top30_bot70 
                    -n1 top30 
                    -n2 bottom70 
                    -c 18
"""

set -euo pipefail

# ============================================
# 使用说明
# ============================================
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

TOBIAS BINDetect comparison between two peak sets

Required Options:
  -1, --peaks1 FILE         First peak set (e.g., Top30 regions)
  -2, --peaks2 FILE         Second peak set (e.g., Bottom70 regions)
  -f, --footprint FILE      Footprint BigWig file
  -m, --motifs FILE         Motif database (JASPAR format)
  -g, --genome FILE         Reference genome FASTA

Optional:
  -o, --outdir DIR          Output directory (default: current directory)
  -p, --prefix PREFIX       Output prefix (default: comparison)
  -n1, --name1 NAME         Name for first condition (default: Condition1)
  -n2, --name2 NAME         Name for second condition (default: Condition2)
  -c, --cores INT           Number of CPU cores (default: 18)
  -h, --help                Show this help message

Example:
  $0 -1 top30.bed -2 bottom70.bed \\
     -f footprint.bw -m motifs.jaspar -g hg19.fa \\
     -o results -p ConEnh -n1 Top30 -n2 Bottom70 -c 18

EOF
    exit 1
}

# ============================================
# 默认参数
# ============================================
OUTDIR="."
PREFIX="comparison"
NAME1="Condition1"
NAME2="Condition2"
CORES=18

PEAKS1=""
PEAKS2=""
FOOTPRINT=""
MOTIFS=""
GENOME=""

# ============================================
# 解析命令行参数
# ============================================
while [[ $# -gt 0 ]]; do
    case $1 in
        -1|--peaks1)
            PEAKS1="$2"
            shift 2
            ;;
        -2|--peaks2)
            PEAKS2="$2"
            shift 2
            ;;
        -f|--footprint)
            FOOTPRINT="$2"
            shift 2
            ;;
        -m|--motifs)
            MOTIFS="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -n1|--name1)
            NAME1="$2"
            shift 2
            ;;
        -n2|--name2)
            NAME2="$2"
            shift 2
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            ;;
    esac
done

# ============================================
# 检查必需参数
# ============================================
if [ -z "$PEAKS1" ] || [ -z "$PEAKS2" ] || [ -z "$FOOTPRINT" ] || \
   [ -z "$MOTIFS" ] || [ -z "$GENOME" ]; then
    echo "ERROR: Missing required arguments!"
    echo ""
    usage
fi

# ============================================
# 创建输出目录
# ============================================
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ============================================
# 显示配置信息
# ============================================
echo "========================================================================"
echo "TOBIAS BINDetect: ${NAME1} vs ${NAME2} Comparison"
echo "========================================================================"
echo ""
echo "Configuration:"
echo "  Peak Set 1:    $PEAKS1"
echo "  Peak Set 2:    $PEAKS2"
echo "  Footprint:     $FOOTPRINT"
echo "  Motifs:        $MOTIFS"
echo "  Genome:        $GENOME"
echo "  Output Dir:    $OUTDIR"
echo "  Output Prefix: $PREFIX"
echo "  Condition 1:   $NAME1"
echo "  Condition 2:   $NAME2"
echo "  CPU Cores:     $CORES"
echo ""

# ============================================
# 检查文件是否存在
# ============================================
echo "Checking input files..."
for file in "$PEAKS1" "$PEAKS2" "$FOOTPRINT" "$MOTIFS" "$GENOME"; do
    if [ ! -f "$file" ]; then
        echo "  ✗ ERROR: File not found: $file"
        exit 1
    fi
    echo "  ✓ $(basename "$file")"
done

echo ""
echo "Peak statistics:"
echo "  ${NAME1}: $(wc -l < "$PEAKS1") regions"
echo "  ${NAME2}: $(wc -l < "$PEAKS2") regions"
echo ""

# ============================================
# 创建临时 header 文件
# ============================================
TEMP_HEADER="${PREFIX}_temp_header.txt"
echo -e "chr\tstart\tend" > "$TEMP_HEADER"

# ============================================
# 定义输出目录
# ============================================
OUTDIR1="${PREFIX}_BINDetect_${NAME1}"
OUTDIR2="${PREFIX}_BINDetect_${NAME2}"

# ============================================
# 运行 BINDetect for Condition 1
# ============================================
echo "========================================================================"
echo "Step 1/3: Running BINDetect for ${NAME1} regions..."
echo "========================================================================"

TOBIAS BINDetect \
    --motifs "$MOTIFS" \
    --signals "$FOOTPRINT" \
    --genome "$GENOME" \
    --peaks "$PEAKS1" \
    --peak_header "$TEMP_HEADER" \
    --outdir "$OUTDIR1" \
    --cond_names "$NAME1" \
    --cores "$CORES"

echo ""
echo "✓ ${NAME1} analysis complete"
echo ""

# ============================================
# 运行 BINDetect for Condition 2
# ============================================
echo "========================================================================"
echo "Step 2/3: Running BINDetect for ${NAME2} regions..."
echo "========================================================================"

TOBIAS BINDetect \
    --motifs "$MOTIFS" \
    --signals "$FOOTPRINT" \
    --genome "$GENOME" \
    --peaks "$PEAKS2" \
    --peak_header "$TEMP_HEADER" \
    --outdir "$OUTDIR2" \
    --cond_names "$NAME2" \
    --cores "$CORES"

echo ""
echo "✓ ${NAME2} analysis complete"
echo ""

# ============================================
# 比较结果
# ============================================
echo "========================================================================"
echo "Step 3/3: Comparing results..."
echo "========================================================================"

python3 << PYEOF
import pandas as pd
import numpy as np
import os
import sys

# 配置
prefix = "${PREFIX}"
name1 = "${NAME1}"
name2 = "${NAME2}"
outdir1 = "${OUTDIR1}"
outdir2 = "${OUTDIR2}"

# 输出文件
comparison_file = f"{prefix}_{name1}_vs_{name2}_comparison.txt"
summary_file = f"{prefix}_{name1}_vs_{name2}_summary.txt"

# 读取结果
print("Loading results...")
file1 = os.path.join(outdir1, "bindetect_results.txt")
file2 = os.path.join(outdir2, "bindetect_results.txt")

if not os.path.exists(file1):
    print(f"ERROR: Result file not found: {file1}")
    sys.exit(1)

if not os.path.exists(file2):
    print(f"ERROR: Result file not found: {file2}")
    sys.exit(1)

data1 = pd.read_csv(file1, sep="\t")
data2 = pd.read_csv(file2, sep="\t")

print(f"  {name1}: {len(data1)} TFs")
print(f"  {name2}: {len(data2)} TFs")

# 合并数据
print("\nMerging data...")
merged = data1.merge(
    data2, 
    on=['name', 'motif_id'], 
    suffixes=(f'_{name1}', f'_{name2}'),
    how='inner'
)

print(f"  Merged: {len(merged)} TFs")

# 查找分数列
score_cols = [c for c in merged.columns if 'mean_score' in c.lower()]
bound_cols = [c for c in merged.columns if 'bound' in c.lower() and 'unbound' not in c.lower()]

print(f"\nFound columns:")
print(f"  Score columns: {score_cols}")
print(f"  Bound columns: {bound_cols}")

if len(score_cols) >= 2:
    # 计算差异
    merged[f'score_{name1}'] = merged[score_cols[0]]
    merged[f'score_{name2}'] = merged[score_cols[1]]
    merged['score_diff'] = merged[f'score_{name1}'] - merged[f'score_{name2}']
    merged['score_log2fc'] = np.log2(
        (merged[f'score_{name1}'] + 0.01) / (merged[f'score_{name2}'] + 0.01)
    )
    
    # 添加 bound 信息
    if len(bound_cols) >= 2:
        merged[f'bound_{name1}'] = merged[bound_cols[0]]
        merged[f'bound_{name2}'] = merged[bound_cols[1]]
        merged['bound_diff'] = merged[f'bound_{name1}'] - merged[f'bound_{name2}']
    
    # 准备输出列
    output_cols = [
        'name', 'motif_id', 
        f'score_{name1}', f'score_{name2}', 
        'score_diff', 'score_log2fc'
    ]
    
    if 'bound_diff' in merged.columns:
        output_cols.extend([
            f'bound_{name1}', f'bound_{name2}', 'bound_diff'
        ])
    
    # 排序并保存
    merged_sorted = merged[output_cols].sort_values('score_diff', ascending=False)
    merged_sorted.to_csv(comparison_file, sep="\t", index=False)
    
    # 显示结果
    print("\n" + "="*80)
    print(f"TOP 10 TFs ENRICHED IN {name1.upper()}")
    print("="*80)
    top10 = merged_sorted.head(10)
    print(top10.to_string(index=False))
    
    print("\n" + "="*80)
    print(f"TOP 10 TFs ENRICHED IN {name2.upper()}")
    print("="*80)
    bottom10 = merged_sorted.tail(10)
    print(bottom10.to_string(index=False))
    
    # 统计信息
    n_enriched_1 = (merged['score_diff'] > 0).sum()
    n_enriched_2 = (merged['score_diff'] < 0).sum()
    mean_diff = merged['score_diff'].mean()
    median_diff = merged['score_diff'].median()
    
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"Total TFs analyzed: {len(merged)}")
    print(f"TFs enriched in {name1} (score_diff > 0): {n_enriched_1} ({n_enriched_1/len(merged)*100:.1f}%)")
    print(f"TFs enriched in {name2} (score_diff < 0): {n_enriched_2} ({n_enriched_2/len(merged)*100:.1f}%)")
    print(f"Mean score difference: {mean_diff:.4f}")
    print(f"Median score difference: {median_diff:.4f}")
    
    # 保存摘要
    with open(summary_file, "w") as f:
        f.write(f"=== {name1} vs {name2} Comparison Summary ===\n\n")
        f.write(f"Analysis Date: $(date)\n")
        f.write(f"Prefix: {prefix}\n\n")
        
        f.write("Input Files:\n")
        f.write(f"  {name1} peaks: ${PEAKS1}\n")
        f.write(f"  {name2} peaks: ${PEAKS2}\n")
        f.write(f"  Footprint: ${FOOTPRINT}\n")
        f.write(f"  Motifs: ${MOTIFS}\n\n")
        
        f.write("Statistics:\n")
        f.write(f"  Total TFs: {len(merged)}\n")
        f.write(f"  Enriched in {name1}: {n_enriched_1} ({n_enriched_1/len(merged)*100:.1f}%)\n")
        f.write(f"  Enriched in {name2}: {n_enriched_2} ({n_enriched_2/len(merged)*100:.1f}%)\n")
        f.write(f"  Mean score difference: {mean_diff:.4f}\n")
        f.write(f"  Median score difference: {median_diff:.4f}\n\n")
        
        f.write(f"Top 10 TFs in {name1}:\n")
        f.write(top10[['name', 'score_diff', 'score_log2fc']].to_string(index=False))
        f.write(f"\n\nTop 10 TFs in {name2}:\n")
        f.write(bottom10[['name', 'score_diff', 'score_log2fc']].to_string(index=False))
    
    print("\n✓ Results saved to:")
    print(f"  - {comparison_file} (full results)")
    print(f"  - {summary_file} (summary)")
    
else:
    print("ERROR: Could not find score columns!")
    print("Available columns:", merged.columns.tolist())
    sys.exit(1)

PYEOF

# ============================================
# 清理临时文件
# ============================================
rm -f "$TEMP_HEADER"

# ============================================
# 最终总结
# ============================================
echo ""
echo "========================================================================"
echo "Analysis Complete!"
echo "========================================================================"
echo "Output directory: $(pwd)"
echo ""
echo "BINDetect results:"
echo "  - ${OUTDIR1}/"
echo "  - ${OUTDIR2}/"
echo ""
echo "Comparison results:"
echo "  - ${PREFIX}_${NAME1}_vs_${NAME2}_comparison.txt"
echo "  - ${PREFIX}_${NAME1}_vs_${NAME2}_summary.txt"
echo ""
echo "========================================================================"
