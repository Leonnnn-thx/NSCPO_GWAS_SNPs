#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOBIAS Footprint Comprehensive Analysis Pipeline
Integration: BigWig to BED -> Score Statistics -> P75 Threshold Filtering -> Region Merging -> Width Analysis

Author: Integrated Pipeline
Date: 2025-11-24
Description: 
    1. BigWig to BED conversion (bw2bed)
    2. Analyze score distribution, calculate P75 threshold
    3. Filter BED file using P75 threshold
    4. Merge adjacent regions using bedtools merge
    5. Statistical analysis of merged region width distribution

Usage:
    python footprint_pipeline.py -i input.bw -o output_prefix
    
Parameters:
    -i, --input: Input BigWig file
    -o, --output: Output file prefix
    -m, --merge_distance: Merge distance for bedtools merge (default: 0 bp)
    --score_bins: Number of bins for score analysis (default: 50)
    --width_bins: Number of bins for width analysis (default: 20)
    --keep_intermediate: Keep intermediate files
    --percentile: Percentile threshold to use (default: 75)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
import argparse
import os
import sys
import subprocess
from datetime import datetime

# ============================================================================
# Global Settings
# ============================================================================
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
sns.set_palette("husl")

# ============================================================================
# Command Line Argument Parsing
# ============================================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description='TOBIAS Footprint Comprehensive Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
    # Basic usage
    python %(prog)s -i footprints.bw -o output
    
    # Specify merge distance
    python %(prog)s -i footprints.bw -o output -m 10
    
    # Use P90 threshold
    python %(prog)s -i footprints.bw -o output --percentile 90
    
    # Keep intermediate files
    python %(prog)s -i footprints.bw -o output --keep_intermediate
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='Input BigWig file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file prefix')
    
    # Optional arguments
    parser.add_argument('-m', '--merge_distance', type=int, default=0,
                       help='Merge distance for bedtools merge (default: 0 bp)')
    parser.add_argument('--score_bins', type=int, default=50,
                       help='Number of bins for score histogram (default: 50)')
    parser.add_argument('--width_bins', type=int, default=20,
                       help='Number of bins for width histogram (default: 20)')
    parser.add_argument('--keep_intermediate', action='store_true',
                       help='Keep intermediate files')
    parser.add_argument('--percentile', type=float, default=75,
                       help='Percentile threshold to use (default: 75)')
    
    return parser.parse_args()

# ============================================================================
# Step 1: BigWig to BED Conversion
# ============================================================================
def bigwig_to_bed(bw_file, bed_file, threshold=0):
    """Convert BigWig file to BED format"""
    print("\n" + "="*80)
    print("[Step 1] BigWig to BED Conversion")
    print("="*80)
    print(f"Input file: {bw_file}")
    print(f"Output file: {bed_file}")
    print(f"Threshold: >= {threshold}")
    print()
    
    try:
        import pyBigWig
    except ImportError:
        print("❌ Error: pyBigWig is required")
        print("   Please run: pip install pyBigWig")
        sys.exit(1)
    
    # Open BigWig file
    bw = pyBigWig.open(bw_file)
    
    total_regions = 0
    filtered_regions = 0
    
    # Output file
    with open(bed_file, 'w') as out:
        # Iterate through all chromosomes
        for chrom in bw.chroms():
            print(f"Processing {chrom}...", end=" ")
            
            # Get all intervals for this chromosome
            intervals = bw.intervals(chrom)
            
            if intervals is None:
                print("No signal")
                continue
            
            # Filter and write
            chr_count = 0
            for start, end, value in intervals:
                total_regions += 1
                
                if value >= threshold:
                    out.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\n")
                    filtered_regions += 1
                    chr_count += 1
            
            print(f"✓ {chr_count:,} regions")
    
    bw.close()
    
    print()
    print(f"{'Conversion Complete':^80}")
    print(f"{'-'*80}")
    print(f"Total regions: {total_regions:,}")
    print(f"After filtering: {filtered_regions:,}")
    print(f"Retention rate: {filtered_regions/total_regions*100:.1f}%")
    print(f"{'-'*80}")
    
    return filtered_regions

# ============================================================================
# Step 2: Score Statistical Analysis
# ============================================================================
def analyze_score_distribution(bed_file, output_prefix, bins=50, percentile=75):
    """Analyze score distribution and calculate threshold"""
    print("\n" + "="*80)
    print("[Step 2] Score Statistical Analysis")
    print("="*80)
    print(f"Reading file: {bed_file}")
    
    # Read data
    df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'score'])
    print(f"✓ Successfully read {len(df):,} sites")
    
    scores = df['score'].values
    
    # Calculate statistics
    stats_dict = {
        'Total Sites': len(scores),
        'Mean': np.mean(scores),
        'Median': np.median(scores),
        'Std Dev': np.std(scores),
        'Min': np.min(scores),
        'Max': np.max(scores),
        '25th Percentile (Q1)': np.percentile(scores, 25),
        '50th Percentile (Q2)': np.percentile(scores, 50),
        '75th Percentile (Q3)': np.percentile(scores, 75),
        '95th Percentile (P95)': np.percentile(scores, 95),
        'Skewness': stats.skew(scores),
        'Kurtosis': stats.kurtosis(scores),
    }
    
    threshold = np.percentile(scores, percentile)
    
    # Print statistics
    print(f"\n{'Score Statistics Summary':^80}")
    print(f"{'-'*80}")
    for key, value in stats_dict.items():
        if isinstance(value, float):
            print(f"{key:.<50} {value:>25.4f}")
        else:
            print(f"{key:.<50} {value:>25,}")
    print(f"{'-'*80}")
    print(f"\n⭐ P{percentile} Threshold: {threshold:.4f}")
    print(f"   Will be used for filtering, retaining top {100-percentile}% high-score sites")
    
    # Generate analysis plots
    print(f"\nGenerating score analysis plots...")
    plot_score_analysis(df, output_prefix, bins, threshold, percentile)
    
    # Save report
    save_score_report(df, stats_dict, output_prefix, threshold, percentile)
    
    return threshold

def plot_score_analysis(df, output_prefix, bins, threshold, percentile):
    """Generate score analysis plots"""
    scores = df['score'].values
    
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Box Plot
    ax1 = plt.subplot(3, 3, 1)
    box_plot = ax1.boxplot(scores, vert=True, patch_artist=True,
                           boxprops=dict(facecolor='lightblue', alpha=0.7),
                           medianprops=dict(color='red', linewidth=2))
    ax1.axhline(threshold, color='red', linestyle='--', linewidth=2,
               label=f'P{percentile}: {threshold:.2f}')
    ax1.set_ylabel('Footprint Score', fontsize=12, fontweight='bold')
    ax1.set_title('Box Plot', fontsize=14, fontweight='bold', pad=15)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Violin Plot
    ax2 = plt.subplot(3, 3, 2)
    parts = ax2.violinplot([scores], vert=True, showmeans=True,
                           showmedians=True, showextrema=True)
    for pc in parts['bodies']:
        pc.set_facecolor('lightcoral')
        pc.set_alpha(0.7)
    ax2.axhline(threshold, color='red', linestyle='--', linewidth=2)
    ax2.set_ylabel('Footprint Score', fontsize=12, fontweight='bold')
    ax2.set_title('Violin Plot', fontsize=14, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.3)
    
    # 3. Histogram - Frequency
    ax3 = plt.subplot(3, 3, 3)
    n, bins_arr, patches = ax3.hist(scores, bins=bins, color='steelblue',
                                     alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.axvline(threshold, color='red', linestyle='--', linewidth=2,
               label=f'P{percentile} Threshold: {threshold:.2f}')
    ax3.axvline(np.mean(scores), color='orange', linestyle='--', linewidth=2,
               label=f'Mean: {np.mean(scores):.2f}')
    ax3.set_xlabel('Footprint Score', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax3.set_title('Histogram - Frequency', fontsize=14, fontweight='bold', pad=15)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Histogram - Density
    ax4 = plt.subplot(3, 3, 4)
    ax4.hist(scores, bins=bins, color='coral', alpha=0.7,
             edgecolor='black', linewidth=0.5, density=True)
    from scipy.stats import gaussian_kde
    density = gaussian_kde(scores)
    xs = np.linspace(scores.min(), scores.max(), 200)
    ax4.plot(xs, density(xs), 'r-', linewidth=2, label='KDE')
    ax4.axvline(threshold, color='blue', linestyle='--', linewidth=2,
               label=f'P{percentile}')
    ax4.set_xlabel('Footprint Score', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax4.set_title('Histogram - Density', fontsize=14, fontweight='bold', pad=15)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. Cumulative Distribution
    ax5 = plt.subplot(3, 3, 5)
    sorted_scores = np.sort(scores)
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores) * 100
    ax5.plot(sorted_scores, cumulative, linewidth=2, color='purple')
    ax5.axvline(threshold, color='red', linestyle='--', linewidth=2)
    ax5.axhline(percentile, color='red', linestyle=':', alpha=0.5)
    ax5.set_xlabel('Footprint Score', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Cumulative Frequency (%)', fontsize=12, fontweight='bold')
    ax5.set_title('Cumulative Distribution (CDF)', fontsize=14, fontweight='bold', pad=15)
    ax5.grid(True, alpha=0.3)
    
    # 6. Log-scale Histogram
    ax6 = plt.subplot(3, 3, 6)
    ax6.hist(scores, bins=bins, color='teal', alpha=0.7,
             edgecolor='black', linewidth=0.5)
    ax6.axvline(threshold, color='red', linestyle='--', linewidth=2)
    ax6.set_xlabel('Footprint Score', fontsize=12, fontweight='bold')
    ax6.set_ylabel('Frequency (log scale)', fontsize=12, fontweight='bold')
    ax6.set_yscale('log')
    ax6.set_title('Histogram - Log Scale', fontsize=14, fontweight='bold', pad=15)
    ax6.grid(True, alpha=0.3, which='both')
    
    # 7. Q-Q Plot
    ax7 = plt.subplot(3, 3, 7)
    stats.probplot(scores, dist="norm", plot=ax7)
    ax7.set_title('Q-Q Plot (Normality Test)', fontsize=14, fontweight='bold', pad=15)
    ax7.grid(True, alpha=0.3)
    
    # 8. ECDF
    ax8 = plt.subplot(3, 3, 8)
    sorted_scores = np.sort(scores)
    y = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    ax8.step(sorted_scores, y, where='post', linewidth=2, color='darkgreen')
    ax8.axvline(threshold, color='red', linestyle='--', linewidth=2)
    ax8.set_xlabel('Footprint Score', fontsize=12, fontweight='bold')
    ax8.set_ylabel('ECDF', fontsize=12, fontweight='bold')
    ax8.set_title('Empirical CDF', fontsize=14, fontweight='bold', pad=15)
    ax8.grid(True, alpha=0.3)
    
    # 9. Statistics Text
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    stats_text = "Statistics Summary\n" + "="*40 + "\n"
    stats_text += f"Total Sites: {len(scores):,}\n"
    stats_text += f"Mean: {np.mean(scores):.4f}\n"
    stats_text += f"Median: {np.median(scores):.4f}\n"
    stats_text += f"Std Dev: {np.std(scores):.4f}\n"
    stats_text += f"Min: {np.min(scores):.4f}\n"
    stats_text += f"Max: {np.max(scores):.4f}\n"
    stats_text += f"\nPercentiles:\n"
    stats_text += f"P25: {np.percentile(scores, 25):.4f}\n"
    stats_text += f"P50: {np.percentile(scores, 50):.4f}\n"
    stats_text += f"P75: {np.percentile(scores, 75):.4f}\n"
    stats_text += f"P95: {np.percentile(scores, 95):.4f}\n"
    stats_text += f"\n⭐ P{percentile} Threshold: {threshold:.4f}\n"
    stats_text += f"\nAfter Filtering:\n"
    stats_text += f"{(scores >= threshold).sum():,} sites\n"
    stats_text += f"Retention: {(scores >= threshold).sum()/len(scores)*100:.1f}%"
    
    ax9.text(0.1, 0.9, stats_text, transform=ax9.transAxes,
            fontsize=11, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    output_file = f'{output_prefix}_step2_score_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Score analysis plot saved: {output_file}")
    plt.close()

def save_score_report(df, stats_dict, output_prefix, threshold, percentile):
    """Save score statistics report"""
    report_file = f'{output_prefix}_step2_score_report.txt'
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("Footprint Score Statistical Analysis Report\n".center(80))
        f.write("="*80 + "\n\n")
        
        f.write(f"Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("1. Basic Statistics\n")
        f.write("-"*80 + "\n")
        for key, value in stats_dict.items():
            if isinstance(value, float):
                f.write(f"{key:.<50} {value:>25.6f}\n")
            else:
                f.write(f"{key:.<50} {value:>25,}\n")
        
        f.write("\n2. Threshold Information\n")
        f.write("-"*80 + "\n")
        f.write(f"Percentile Used: P{percentile}\n")
        f.write(f"Threshold: {threshold:.6f}\n")
        f.write(f"Sites Retained After Filtering: {(df['score'] >= threshold).sum():,}\n")
        f.write(f"Retention Rate: {(df['score'] >= threshold).sum()/len(df)*100:.2f}%\n")
        
        f.write("\n" + "="*80 + "\n")
    
    print(f"✓ Score statistics report saved: {report_file}")

# ============================================================================
# Step 3: Filter BED File by Score
# ============================================================================
def filter_bed_by_score(input_bed, output_bed, threshold):
    """Filter BED file by score threshold"""
    print("\n" + "="*80)
    print("[Step 3] Filter BED File by Threshold")
    print("="*80)
    print(f"Input file: {input_bed}")
    print(f"Output file: {output_bed}")
    print(f"Score threshold: >= {threshold:.4f}")
    print()
    
    # Read and filter
    df = pd.read_csv(input_bed, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'score'])
    
    before_count = len(df)
    df_filtered = df[df['score'] >= threshold]
    after_count = len(df_filtered)
    
    # Save filtered file
    df_filtered.to_csv(output_bed, sep='\t', index=False, header=False)
    
    print(f"{'Filtering Complete':^80}")
    print(f"{'-'*80}")
    print(f"Before filtering: {before_count:,} regions")
    print(f"After filtering: {after_count:,} regions")
    print(f"Retention rate: {after_count/before_count*100:.1f}%")
    print(f"{'-'*80}")
    
    return after_count

# ============================================================================
# Step 4: bedtools merge
# ============================================================================
def merge_bed_regions(input_bed, output_bed, merge_distance=0):
    """Merge adjacent regions using bedtools merge"""
    print("\n" + "="*80)
    print("[Step 4] Merge Adjacent Regions")
    print("="*80)
    print(f"Input file: {input_bed}")
    print(f"Output file: {output_bed}")
    print(f"Merge distance: {merge_distance} bp")
    print()
    
    # Check bedtools
    try:
        subprocess.run(['bedtools', '--version'],
                      capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ Error: bedtools not found")
        print("   Please install: conda install -c bioconda bedtools")
        sys.exit(1)
    
    # Sort
    sorted_bed = input_bed.replace('.bed', '_sorted.bed')
    print("Sorting BED file...")
    sort_cmd = f"sort -k1,1 -k2,2n {input_bed} > {sorted_bed}"
    subprocess.run(sort_cmd, shell=True, check=True)
    print(f"✓ Sorting complete")
    
    # Merge
    print(f"Merging regions...")
    merge_cmd = f"bedtools merge -i {sorted_bed} -d {merge_distance} > {output_bed}"
    subprocess.run(merge_cmd, shell=True, check=True)
    
    # Statistics
    with open(input_bed) as f:
        before_count = sum(1 for _ in f)
    
    with open(output_bed) as f:
        after_count = sum(1 for _ in f)
    
    print()
    print(f"{'Merging Complete':^80}")
    print(f"{'-'*80}")
    print(f"Before merging: {before_count:,} regions")
    print(f"After merging: {after_count:,} regions")
    print(f"Merge rate: {(1 - after_count/before_count)*100:.1f}%")
    print(f"{'-'*80}")
    
    # Clean up temporary files
    if os.path.exists(sorted_bed):
        os.remove(sorted_bed)
    
    return after_count

# ============================================================================
# Step 5: Width Analysis
# ============================================================================
def analyze_width_distribution(bed_file, output_prefix, bins=20):
    """Analyze region width distribution"""
    print("\n" + "="*80)
    print("[Step 5] Region Width Analysis")
    print("="*80)
    print(f"Reading file: {bed_file}")
    
    # Read data
    df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chr', 'start', 'end'])
    df['width'] = df['end'] - df['start']
    
    print(f"✓ Successfully read {len(df):,} regions")
    
    widths = df['width'].values
    
    # Calculate statistics
    stats_dict = {
        'Total Regions': len(widths),
        'Total Length (bp)': np.sum(widths),
        'Mean Width': np.mean(widths),
        'Median': np.median(widths),
        'Std Dev': np.std(widths),
        'Min': np.min(widths),
        'Max': np.max(widths),
        '25th Percentile': np.percentile(widths, 25),
        '75th Percentile': np.percentile(widths, 75),
        'IQR': np.percentile(widths, 75) - np.percentile(widths, 25),
    }
    
    # Print statistics
    print(f"\n{'Width Statistics Summary':^80}")
    print(f"{'-'*80}")
    for key, value in stats_dict.items():
        if isinstance(value, float):
            print(f"{key:.<50} {value:>25,.2f} bp")
        else:
            print(f"{key:.<50} {value:>25,}")
    print(f"{'-'*80}")
    
    # Generate analysis plots
    print(f"\nGenerating width analysis plots...")
    plot_width_analysis(df, output_prefix, bins)
    
    # Save report
    save_width_report(df, stats_dict, output_prefix, bins)
    
    return stats_dict

def plot_width_analysis(df, output_prefix, bins):
    """Generate width analysis plots"""
    widths = df['width'].values
    
    fig = plt.figure(figsize=(18, 12))
    
    # 1. Box Plot
    ax1 = plt.subplot(2, 3, 1)
    box_plot = ax1.boxplot(widths, vert=True, patch_artist=True,
                           boxprops=dict(facecolor='lightblue', alpha=0.7),
                           medianprops=dict(color='red', linewidth=2))
    ax1.set_ylabel('Region Width (bp)', fontsize=12, fontweight='bold')
    ax1.set_title('Box Plot', fontsize=14, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3)
    
    stats_text = f"Median: {np.median(widths):.0f} bp\n"
    stats_text += f"Q1: {np.percentile(widths, 25):.0f} bp\n"
    stats_text += f"Q3: {np.percentile(widths, 75):.0f} bp"
    ax1.text(0.98, 0.97, stats_text, transform=ax1.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            fontsize=10)
    
    # 2. Violin Plot
    ax2 = plt.subplot(2, 3, 2)
    parts = ax2.violinplot([widths], vert=True, showmeans=True,
                           showmedians=True, showextrema=True)
    for pc in parts['bodies']:
        pc.set_facecolor('lightcoral')
        pc.set_alpha(0.7)
    ax2.set_ylabel('Region Width (bp)', fontsize=12, fontweight='bold')
    ax2.set_title('Violin Plot', fontsize=14, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.3)
    
    # 3. Histogram
    ax3 = plt.subplot(2, 3, 3)
    n, bins_arr, patches = ax3.hist(widths, bins=50, color='steelblue',
                                     alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.axvline(np.mean(widths), color='red', linestyle='--', linewidth=2,
               label=f'Mean: {np.mean(widths):.0f} bp')
    ax3.axvline(np.median(widths), color='green', linestyle='--', linewidth=2,
               label=f'Median: {np.median(widths):.0f} bp')
    ax3.set_xlabel('Region Width (bp)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax3.set_title('Histogram', fontsize=14, fontweight='bold', pad=15)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Cumulative Distribution
    ax4 = plt.subplot(2, 3, 4)
    sorted_widths = np.sort(widths)
    cumulative = np.arange(1, len(sorted_widths) + 1) / len(sorted_widths) * 100
    ax4.plot(sorted_widths, cumulative, linewidth=2, color='purple')
    for p in [25, 50, 75]:
        val = np.percentile(widths, p)
        ax4.axvline(val, color='red', linestyle=':', alpha=0.5)
        ax4.axhline(p, color='red', linestyle=':', alpha=0.5)
    ax4.set_xlabel('Region Width (bp)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Cumulative Frequency (%)', fontsize=12, fontweight='bold')
    ax4.set_title('Cumulative Distribution', fontsize=14, fontweight='bold', pad=15)
    ax4.grid(True, alpha=0.3)
    
    # 5. Width Range Bar Chart
    ax5 = plt.subplot(2, 3, 5)
    bins_group = [0, 100, 500, 1000, 2000, 5000, np.inf]
    labels_group = ['0-100', '100-500', '500-1K', '1K-2K', '2K-5K', '>5K']
    df_copy = df.copy()
    df_copy['width_group'] = pd.cut(df_copy['width'], bins=bins_group, labels=labels_group)
    group_counts = df_copy['width_group'].value_counts().sort_index()
    
    bars = ax5.bar(range(len(group_counts)), group_counts.values,
                   color=sns.color_palette('viridis', len(group_counts)),
                   edgecolor='black', linewidth=1)
    ax5.set_xticks(range(len(group_counts)))
    ax5.set_xticklabels(group_counts.index, rotation=45)
    ax5.set_xlabel('Width Range (bp)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax5.set_title('Width Distribution by Range', fontsize=14, fontweight='bold', pad=15)
    ax5.grid(True, alpha=0.3, axis='y')
    
    for bar, count in zip(bars, group_counts.values):
        height = bar.get_height()
        percentage = count / len(df) * 100
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                f'{percentage:.1f}%',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 6. Statistics Text
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    stats_text = "Width Statistics Summary\n" + "="*40 + "\n"
    stats_text += f"Total Regions: {len(widths):,}\n"
    stats_text += f"Total Length: {np.sum(widths):,.0f} bp\n"
    stats_text += f"Mean Width: {np.mean(widths):.2f} bp\n"
    stats_text += f"Median: {np.median(widths):.2f} bp\n"
    stats_text += f"Std Dev: {np.std(widths):.2f} bp\n"
    stats_text += f"Min: {np.min(widths):.0f} bp\n"
    stats_text += f"Max: {np.max(widths):.0f} bp\n"
    stats_text += f"\nPercentiles:\n"
    stats_text += f"P25: {np.percentile(widths, 25):.2f} bp\n"
    stats_text += f"P50: {np.percentile(widths, 50):.2f} bp\n"
    stats_text += f"P75: {np.percentile(widths, 75):.2f} bp\n"
    stats_text += f"\nWidth Groups:\n"
    for label, count in group_counts.items():
        pct = count / len(widths) * 100
        stats_text += f"{label} bp: {count:,} ({pct:.1f}%)\n"
    
    ax6.text(0.1, 0.9, stats_text, transform=ax6.transAxes,
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    plt.tight_layout()
    output_file = f'{output_prefix}_step5_width_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Width analysis plot saved: {output_file}")
    plt.close()

def save_width_report(df, stats_dict, output_prefix, bins):
    """Save width statistics report"""
    report_file = f'{output_prefix}_step5_width_report.txt'
    
    widths = df['width'].values
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("Region Width Statistical Analysis Report\n".center(80))
        f.write("="*80 + "\n\n")
        
        f.write(f"Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("1. Basic Statistics\n")
        f.write("-"*80 + "\n")
        for key, value in stats_dict.items():
            if isinstance(value, float):
                f.write(f"{key:.<50} {value:>25,.2f} bp\n")
            else:
                f.write(f"{key:.<50} {value:>25,}\n")
        
        f.write("\n2. Width Range Group Statistics\n")
        f.write("-"*80 + "\n")
        bins_group = [0, 100, 500, 1000, 2000, 5000, np.inf]
        labels_group = ['0-100bp', '100-500bp', '500-1Kbp', '1K-2Kbp', '2K-5Kbp', '>5Kbp']
        df_copy = df.copy()
        df_copy['width_group'] = pd.cut(df_copy['width'], bins=bins_group, labels=labels_group)
        group_stats = df_copy['width_group'].value_counts().sort_index()
        
        for label, count in group_stats.items():
            pct = count / len(df) * 100
            f.write(f"{label:.<50} {count:>10,} ({pct:>6.2f}%)\n")
        
        f.write("\n" + "="*80 + "\n")
    
    print(f"✓ Width statistics report saved: {report_file}")
    
    # Save final BED file
    output_bed = f'{output_prefix}_final_with_width.bed'
    df.to_csv(output_bed, sep='\t', index=False, header=False)
    print(f"✓ Final BED file saved: {output_bed}")

# ============================================================================
# Main Pipeline
# ============================================================================
def main():
    """Main function - integrate all steps"""
    
    # Parse arguments
    args = parse_args()
    
    # Print pipeline information
    print("\n" + "="*80)
    print("TOBIAS Footprint Comprehensive Analysis Pipeline".center(80))
    print("="*80)
    print(f"\nPipeline Overview:")
    print(f"  1. BigWig to BED Conversion")
    print(f"  2. Score Statistical Analysis -> Calculate P{args.percentile} Threshold")
    print(f"  3. Filter BED by P{args.percentile} Threshold")
    print(f"  4. Merge Adjacent Regions using bedtools merge (distance {args.merge_distance} bp)")
    print(f"  5. Width Statistical Analysis")
    print(f"\nInput file: {args.input}")
    print(f"Output prefix: {args.output}")
    print(f"Keep intermediate files: {'Yes' if args.keep_intermediate else 'No'}")
    print("="*80)
    
    # Define intermediate file names
    step1_bed = f'{args.output}_step1_raw.bed'
    step3_bed = f'{args.output}_step3_filtered.bed'
    step4_bed = f'{args.output}_step4_merged.bed'
    
    try:
        # ========== Step 1: BigWig to BED ==========
        bigwig_to_bed(args.input, step1_bed, threshold=0)
        
        # ========== Step 2: Score Statistical Analysis ==========
        threshold = analyze_score_distribution(
            step1_bed, args.output, 
            bins=args.score_bins,
            percentile=args.percentile
        )
        
        # ========== Step 3: Filter BED File ==========
        filter_bed_by_score(step1_bed, step3_bed, threshold)
        
        # ========== Step 4: bedtools merge ==========
        merge_bed_regions(step3_bed, step4_bed, args.merge_distance)
        
        # ========== Step 5: Width Statistical Analysis ==========
        analyze_width_distribution(step4_bed, args.output, bins=args.width_bins)
        
        # ========== Clean up intermediate files ==========
        if not args.keep_intermediate:
            print(f"\nCleaning up intermediate files...")
            for f in [step1_bed, step3_bed]:
                if os.path.exists(f):
                    os.remove(f)
                    print(f"✓ Deleted: {f}")
        
        # ========== Generate final summary ==========
        print("\n" + "="*80)
        print("Pipeline Complete!".center(80))
        print("="*80)
        print("\nGenerated Files:")
        print(f"\n[Step 1] BigWig to BED:")
        if args.keep_intermediate:
            print(f"  - {step1_bed}")
        else:
            print(f"  - (cleaned up)")
        print(f"\n[Step 2] Score Analysis:")
        print(f"  - {args.output}_step2_score_analysis.png")
        print(f"  - {args.output}_step2_score_report.txt")
        print(f"\n[Step 3] Filtered BED:")
        if args.keep_intermediate:
            print(f"  - {step3_bed}")
        else:
            print(f"  - (cleaned up)")
        print(f"\n[Step 4] Merged BED:")
        print(f"  - {step4_bed}")
        print(f"\n[Step 5] Width Analysis:")
        print(f"  - {args.output}_step5_width_analysis.png")
        print(f"  - {args.output}_step5_width_report.txt")
        print(f"  - {args.output}_final_with_width.bed")
        
        print("\n" + "="*80)
        print("Key Parameters:")
        print(f"  P{args.percentile} Threshold: {threshold:.4f}")
        print(f"  Merge Distance: {args.merge_distance} bp")
        print("="*80 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
