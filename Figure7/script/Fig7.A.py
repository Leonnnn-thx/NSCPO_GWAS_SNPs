import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import argparse
import sys
import os
from pybedtools import BedTool
import pybedtools

# ============================================
# Command Line Argument Parsing
# ============================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Plot Footprint Score Cumulative Distribution Function (CDF) and extract bottom percentile regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Use default parameters (70th percentile)
  python script.py -i input.bed -o output_prefix
  
  # Custom percentile (80th)
  python script.py -i input.bed -o output_prefix -p 80
  
  # Output PDF only
  python script.py -i input.bed -o output_prefix --pdf-only
  
  # Output PNG only
  python script.py -i input.bed -o output_prefix --png-only
  
  # Extract bottom 30% and merge
  python script.py -i input.bed -o output_prefix -p 70 --extract-bottom --merge
        """
    )
    
    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input BED file path (required)')
    
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output file prefix (required)')
    
    parser.add_argument('-p', '--percentile',
                        type=float,
                        default=70.0,
                        help='Percentile threshold (0-100), default: 70')
    
    parser.add_argument('--pdf-only',
                        action='store_true',
                        help='Output PDF format only')
    
    parser.add_argument('--png-only',
                        action='store_true',
                        help='Output PNG format only')
    
    parser.add_argument('--dpi',
                        type=int,
                        default=300,
                        help='PNG output resolution, default: 300')
    
    parser.add_argument('--color',
                        type=str,
                        default='#E74C3C',
                        help='Curve color (hex code), default: #E74C3C (red)')
    
    parser.add_argument('--title',
                        type=str,
                        default='Cumulative Distribution of Footprint Score',
                        help='Plot title')
    
    parser.add_argument('--extract-bottom',
                        action='store_true',
                        help='Extract bottom percentile regions to BED file')
    
    parser.add_argument('--merge',
                        action='store_true',
                        help='Merge overlapping regions in bottom percentile BED file (requires --extract-bottom)')
    
    parser.add_argument('--merge-distance',
                        type=int,
                        default=0,
                        help='Maximum distance between features to merge, default: 0 (only overlapping)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        parser.error(f"Input file does not exist: {args.input}")
    
    # Validate percentile
    if not 0 < args.percentile < 100:
        parser.error(f"Percentile must be between 0-100, current value: {args.percentile}")
    
    # Validate output format
    if args.pdf_only and args.png_only:
        parser.error("Cannot specify both --pdf-only and --png-only")
    
    # Validate merge option
    if args.merge and not args.extract_bottom:
        parser.error("--merge requires --extract-bottom")
    
    return args

# ============================================
# Data Loading Function
# ============================================

def load_data(file_path):
    """Load BED format data"""
    try:
        data = pd.read_csv(file_path, 
                          sep='\t', 
                          header=None,
                          names=['chr', 'start', 'end', 'footprint_score'])
        
        # Check for missing values
        if data['footprint_score'].isna().any():
            n_missing = data['footprint_score'].isna().sum()
            print(f"⚠️  Warning: Data contains {n_missing} missing values, will be removed")
            data = data.dropna(subset=['footprint_score'])
        
        return data
    
    except Exception as e:
        print(f"❌ Failed to read file: {e}")
        sys.exit(1)

# ============================================
# Statistical Analysis Functions
# ============================================

def calculate_statistics(data, percentile):
    """Calculate statistical information"""
    stats = {
        'count': len(data),
        'mean': data['footprint_score'].mean(),
        'std': data['footprint_score'].std(),
        'min': data['footprint_score'].min(),
        'max': data['footprint_score'].max(),
        'median': data['footprint_score'].median(),
        'percentile_value': np.percentile(data['footprint_score'], percentile)
    }
    return stats

def print_statistics(stats, percentile):
    """Print statistical information"""
    print("="*60)
    print("Data Statistics")
    print("="*60)
    print(f"Total rows:           {stats['count']:,}")
    print(f"Mean:                 {stats['mean']:.4f}")
    print(f"Std Dev:              {stats['std']:.4f}")
    print(f"Min:                  {stats['min']:.4f}")
    print(f"Max:                  {stats['max']:.4f}")
    print(f"Median:               {stats['median']:.4f}")
    print(f"{percentile}% Percentile:    {stats['percentile_value']:.4f}")
    print("="*60 + "\n")

# ============================================
# CDF Calculation Function
# ============================================

def calculate_cdf(data):
    """Calculate cumulative distribution function"""
    sorted_data = data.sort_values('footprint_score').reset_index(drop=True)
    sorted_data['cumulative_freq'] = (np.arange(1, len(sorted_data) + 1) / len(sorted_data)) * 100
    return sorted_data

# ============================================
# Bottom Percentile Extraction Function
# ============================================

def extract_bottom_percentile(data, percentile, output_prefix, merge=False, merge_distance=0):
    """Extract bottom percentile regions and optionally merge them"""
    
    # Calculate threshold
    threshold = np.percentile(data['footprint_score'], percentile)
    
    # Filter bottom percentile
    bottom_data = data[data['footprint_score'] <= threshold].copy()
    bottom_data = bottom_data.sort_values(['chr', 'start'])
    
    print(f"\nExtracting bottom {percentile:.1f}% regions...")
    print(f"Threshold: {threshold:.4f}")
    print(f"Number of regions: {len(bottom_data):,}")
    
    # Save unmerged BED file
    unmerged_bed = f"{output_prefix}_bottom{int(percentile)}.bed"
    bottom_data[['chr', 'start', 'end', 'footprint_score']].to_csv(
        unmerged_bed, 
        sep='\t', 
        header=False, 
        index=False
    )
    print(f"✓ Saved unmerged BED: {unmerged_bed}")
    
    saved_files = [unmerged_bed]
    
    # Merge if requested
    if merge:
        print(f"\nMerging overlapping regions (distance: {merge_distance})...")
        
        try:
            # Create BedTool object
            bed = BedTool(unmerged_bed)
            
            # Sort and merge
            merged_bed = bed.sort().merge(d=merge_distance)
            
            # Save merged BED file
            merged_bed_file = f"{output_prefix}_bottom{int(percentile)}_merged.bed"
            merged_bed.saveas(merged_bed_file)
            
            # Count merged regions
            n_merged = len(merged_bed)
            reduction = ((len(bottom_data) - n_merged) / len(bottom_data)) * 100
            
            print(f"✓ Saved merged BED: {merged_bed_file}")
            print(f"  Original regions: {len(bottom_data):,}")
            print(f"  Merged regions:   {n_merged:,}")
            print(f"  Reduction:        {reduction:.1f}%")
            
            saved_files.append(merged_bed_file)
            
            # Calculate statistics for merged regions
            merged_df = pd.read_csv(merged_bed_file, sep='\t', header=None,
                                   names=['chr', 'start', 'end'])
            merged_df['length'] = merged_df['end'] - merged_df['start']
            
            print(f"\nMerged region statistics:")
            print(f"  Total length:     {merged_df['length'].sum():,} bp")
            print(f"  Mean length:      {merged_df['length'].mean():.1f} bp")
            print(f"  Median length:    {merged_df['length'].median():.1f} bp")
            print(f"  Max length:       {merged_df['length'].max():,} bp")
            
        except Exception as e:
            print(f"❌ Failed to merge regions: {e}")
            print("   Make sure pybedtools is properly installed")
    
    return saved_files, bottom_data

# ============================================
# Plotting Function
# ============================================

def plot_cdf(sorted_data, stats, percentile, color, title, bottom_threshold=None):
    """Plot CDF curve"""
    
    # Set plot style
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    # Create figure
    fig, ax = plt.subplots(figsize=(9, 7))
    
    # Plot CDF curve
    ax.plot(sorted_data['footprint_score'], 
            sorted_data['cumulative_freq'],
            color=color,
            linewidth=2.5,
            zorder=3,
            label='Footprint Score CDF')
    
    # Add percentile horizontal line
    ax.axhline(y=percentile, 
               linestyle='--', 
               color='black', 
               linewidth=1,
               alpha=0.5,
               zorder=2,
               label=f'{percentile}% Percentile')
    
    # Add percentile vertical line
    ax.axvline(x=stats['percentile_value'], 
               linestyle='--', 
               color=color, 
               linewidth=1.2,
               alpha=0.7,
               zorder=2)
    
    # Add bottom percentile shading if threshold provided
    if bottom_threshold is not None:
        ax.axvspan(0, bottom_threshold, 
                   alpha=0.2, 
                   color=color,
                   zorder=1,
                   label=f'Bottom {percentile:.0f}%')
    
    # Add percentile value annotation
    ax.text(stats['percentile_value'], 5, 
            f'{stats["percentile_value"]:.3f}',
            fontsize=11,
            color=color,
            ha='center',
            va='bottom',
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', 
                      facecolor='white', 
                      edgecolor=color, 
                      alpha=0.9,
                      linewidth=1.5))
    
    # Add sample size annotation
    ax.text(0.98, 0.02,
            f'n = {stats["count"]:,}',
            transform=ax.transAxes,
            fontsize=11,
            ha='right',
            va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', 
                      facecolor='white', 
                      edgecolor='gray', 
                      alpha=0.9))
    
    # Set axes
    ax.set_xlabel('Footprint Score', fontsize=14, fontweight='normal')
    ax.set_ylabel('Cumulative Frequency (%)', fontsize=14, fontweight='normal')
    ax.set_title(title, 
                 fontsize=16, 
                 fontweight='bold',
                 pad=15)
    
    # Set x-axis range
    x_max = stats['max']
    ax.set_xlim(0, x_max * 1.02)
    
    # Set y-axis range and ticks
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    
    # Add grid
    ax.grid(True, which='major', color='lightgray', linewidth=0.5, alpha=0.5, zorder=1)
    ax.set_axisbelow(True)
    
    # Set borders
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.2)
    
    # Adjust tick label size
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Add legend
    ax.legend(loc='lower right', 
              fontsize=12, 
              frameon=True, 
              shadow=True,
              fancybox=True,
              framealpha=0.95)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig

# ============================================
# Save Figure Function
# ============================================

def save_figure(fig, output_prefix, pdf_only, png_only, dpi):
    """Save figure as PDF and/or PNG format"""
    
    saved_files = []
    
    # Save PDF
    if not png_only:
        pdf_path = f"{output_prefix}_CDF.pdf"
        fig.savefig(pdf_path, 
                    format='pdf', 
                    dpi=300,
                    bbox_inches='tight',
                    pad_inches=0.3)
        saved_files.append(pdf_path)
        print(f"✓ Saved PDF: {pdf_path}")
    
    # Save PNG
    if not pdf_only:
        png_path = f"{output_prefix}_CDF.png"
        fig.savefig(png_path, 
                    format='png', 
                    dpi=dpi,
                    bbox_inches='tight',
                    pad_inches=0.3)
        saved_files.append(png_path)
        print(f"✓ Saved PNG: {png_path}")
    
    return saved_files

# ============================================
# Main Function
# ============================================

def main():
    """Main program"""
    
    # Parse command line arguments
    args = parse_arguments()
    
    print("\n" + "="*60)
    print("Footprint Score CDF Analysis Tool")
    print("="*60)
    print(f"Input file:       {args.input}")
    print(f"Output prefix:    {args.output}")
    print(f"Percentile:       {args.percentile}%")
    print(f"Curve color:      {args.color}")
    print(f"Plot title:       {args.title}")
    print(f"Extract bottom:   {args.extract_bottom}")
    if args.extract_bottom:
        print(f"Merge regions:    {args.merge}")
        if args.merge:
            print(f"Merge distance:   {args.merge_distance} bp")
    
    if args.pdf_only:
        print(f"Output format:    PDF only")
    elif args.png_only:
        print(f"Output format:    PNG only (DPI: {args.dpi})")
    else:
        print(f"Output format:    PDF + PNG (DPI: {args.dpi})")
    print("="*60 + "\n")
    
    # Load data
    print("Loading data...")
    data = load_data(args.input)
    print(f"✓ Successfully loaded {len(data):,} rows\n")
    
    # Calculate statistics
    stats = calculate_statistics(data, args.percentile)
    print_statistics(stats, args.percentile)
    
    # Calculate CDF
    print("Calculating cumulative distribution...")
    sorted_data = calculate_cdf(data)
    print("✓ Done\n")
    
    # Extract bottom percentile if requested
    all_saved_files = []
    bottom_threshold = None
    if args.extract_bottom:
        bed_files, bottom_data = extract_bottom_percentile(
            data, 
            args.percentile, 
            args.output, 
            merge=args.merge,
            merge_distance=args.merge_distance
        )
        all_saved_files.extend(bed_files)
        bottom_threshold = np.percentile(data['footprint_score'], args.percentile)
    
    # Plot figure
    print("\nGenerating plot...")
    fig = plot_cdf(sorted_data, stats, args.percentile, args.color, args.title, bottom_threshold)
    print("✓ Done\n")
    
    # Save figure
    print("Saving plot...")
    plot_files = save_figure(fig, args.output, args.pdf_only, args.png_only, args.dpi)
    all_saved_files.extend(plot_files)
    
    # Close figure
    plt.close()
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"Total files saved: {len(all_saved_files)}")
    for file in all_saved_files:
        print(f"  • {file}")
    print("="*60 + "\n")

# ============================================
# Program Entry Point
# ============================================

if __name__ == "__main__":
    main()
