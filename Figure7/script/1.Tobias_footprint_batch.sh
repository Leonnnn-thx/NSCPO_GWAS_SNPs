#!/bin/bash
#############################################################
# TOBIAS Pipeline for ATAC-seq Footprinting Analysis
# Author: Thx
# Date: 2025-11-22
# Description: Automated pipeline for ATACorrect, FootprintScores, and BINDetect
# Usage: ./tobias_pipeline.sh --config config.txt
#        or with command line arguments
#./tobias_pipeline.sh \
#   --genome /path/to/genome.fa \
#    --peaks /path/to/peaks.bed \
#    --motifs /path/to/motifs.jaspar \
#    --samples /path/to/samples.txt \
#    --blacklist /path/to/blacklist.bed \
#    --peak-header /path/to/header.txt \
#    --cores 32 \
#    --outdir /path/to/output
#############################################################

set -euo pipefail  # Exit on error

#############################################################
# Default Configuration
#############################################################

# Reference files (shared across samples)
GENOME="/home/thx/reference_genome/hg19/hg19.fa"
BLACKLIST="/home/thx/reference_genome/hg19/hg19_blacklist_ENCFF001TDO.bed"
PEAKS="/home/thx/NSCPO_Variant/HEPM_ATAC/2.Process_data/footprint/805Enhancer.bed"
PEAK_HEADER="/home/thx/NSCPO_Variant/HEPM_ATAC/2.Process_data/footprint/peaks_header.txt"
MOTIFS="/home/thx/NSCPO_Variant/HEPM_ATAC/2.Process_data/footprint/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar/JASPAR2024_CORE_vertebrates_merged.jaspar"

# Computing resources
CORES=18

# Output base directory
OUTPUT_BASE="/home/thx/NSCPO_Variant/HEPM_ATAC/2.Process_data/footprint"

# Sample list file
SAMPLE_LIST=""

# Individual samples array
declare -a SAMPLES

#############################################################
# Help Message
#############################################################

show_help() {
    cat << EOF
TOBIAS Pipeline for ATAC-seq Footprinting Analysis

Usage: 
    $0 [OPTIONS]

Required Arguments:
    --genome PATH           Path to reference genome FASTA file
    --peaks PATH            Path to peak regions BED file
    --motifs PATH           Path to motif database (JASPAR format)
    --samples PATH          Path to sample list file (TSV format)
                           Format: SAMPLE_NAME<TAB>BAM_PATH

Optional Arguments:
    --blacklist PATH        Path to blacklist regions BED file
    --peak-header PATH      Path to peak header file
    --outdir PATH           Output base directory (default: current directory)
    --cores INT             Number of CPU cores to use (default: 18)
    --config PATH           Path to configuration file
    --skip-atacorrect       Skip ATACorrect step (use existing corrected bigWig)
    --skip-footprint        Skip FootprintScores step
    --skip-bindetect        Skip BINDetect step
    -h, --help              Show this help message

Configuration File Format (config.txt):
    GENOME=/path/to/genome.fa
    BLACKLIST=/path/to/blacklist.bed
    PEAKS=/path/to/peaks.bed
    PEAK_HEADER=/path/to/header.txt
    MOTIFS=/path/to/motifs.jaspar
    CORES=18
    OUTPUT_BASE=/path/to/output
    SAMPLE_LIST=/path/to/samples.txt

Sample List File Format (samples.txt):
    HEPM_ATAC_1<TAB>/path/to/HEPM_ATAC_1_sorted.bam
    HEPM_ATAC_2<TAB>/path/to/HEPM_ATAC_2_sorted.bam

Examples:
    # Using configuration file
    $0 --config config.txt

    # Using command line arguments
    $0 --genome hg19.fa --peaks peaks.bed --motifs motifs.jaspar --samples samples.txt

    # Override config file with command line
    $0 --config config.txt --cores 32 --outdir /new/output/path

EOF
    exit 0
}

#############################################################
# Parse Command Line Arguments
#############################################################

SKIP_ATACORRECT=false
SKIP_FOOTPRINT=false
SKIP_BINDETECT=false

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --genome)
                GENOME="$2"
                shift 2
                ;;
            --blacklist)
                BLACKLIST="$2"
                shift 2
                ;;
            --peaks)
                PEAKS="$2"
                shift 2
                ;;
            --peak-header)
                PEAK_HEADER="$2"
                shift 2
                ;;
            --motifs)
                MOTIFS="$2"
                shift 2
                ;;
            --cores)
                CORES="$2"
                shift 2
                ;;
            --outdir)
                OUTPUT_BASE="$2"
                shift 2
                ;;
            --samples)
                SAMPLE_LIST="$2"
                shift 2
                ;;
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --skip-atacorrect)
                SKIP_ATACORRECT=true
                shift
                ;;
            --skip-footprint)
                SKIP_FOOTPRINT=true
                shift
                ;;
            --skip-bindetect)
                SKIP_BINDETECT=true
                shift
                ;;
            -h|--help)
                show_help
                ;;
            *)
                echo "ERROR: Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
}

#############################################################
# Load Configuration File
#############################################################

load_config() {
    if [ -n "$CONFIG_FILE" ]; then
        if [ ! -f "$CONFIG_FILE" ]; then
            echo "ERROR: Configuration file not found: $CONFIG_FILE"
            exit 1
        fi
        
        echo "Loading configuration from: $CONFIG_FILE"
        
        while IFS='=' read -r key value; do
            # Skip comments and empty lines
            [[ "$key" =~ ^#.*$ ]] && continue
            [[ -z "$key" ]] && continue
            
            # Remove leading/trailing whitespace
            key=$(echo "$key" | xargs)
            value=$(echo "$value" | xargs)
            
            case "$key" in
                GENOME) GENOME="$value" ;;
                BLACKLIST) BLACKLIST="$value" ;;
                PEAKS) PEAKS="$value" ;;
                PEAK_HEADER) PEAK_HEADER="$value" ;;
                MOTIFS) MOTIFS="$value" ;;
                CORES) CORES="$value" ;;
                OUTPUT_BASE) OUTPUT_BASE="$value" ;;
                SAMPLE_LIST) SAMPLE_LIST="$value" ;;
            esac
        done < "$CONFIG_FILE"
    fi
}

#############################################################
# Load Sample List
#############################################################

load_samples() {
    if [ -n "$SAMPLE_LIST" ]; then
        if [ ! -f "$SAMPLE_LIST" ]; then
            echo "ERROR: Sample list file not found: $SAMPLE_LIST"
            exit 1
        fi
        
        echo "Loading samples from: $SAMPLE_LIST"
        
        SAMPLES=()
        while IFS=$'\t' read -r sample_name bam_path; do
            # Skip comments and empty lines
            [[ "$sample_name" =~ ^#.*$ ]] && continue
            [[ -z "$sample_name" ]] && continue
            
            SAMPLES+=("${sample_name}|${bam_path}")
        done < "$SAMPLE_LIST"
        
        if [ ${#SAMPLES[@]} -eq 0 ]; then
            echo "ERROR: No samples found in $SAMPLE_LIST"
            exit 1
        fi
    fi
}

#############################################################
# Functions
#############################################################

print_header() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "============================================================"
    echo ""
}

print_step() {
    echo ""
    echo "------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "------------------------------------------------------------"
}

check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: File not found: $1"
        exit 1
    fi
}

#############################################################
# Validation
#############################################################

validate_config() {
    print_header "Validating Configuration"
    
    local errors=0
    
    # Check required files
    if [ ! -f "$GENOME" ]; then
        echo "✗ ERROR: Genome file not found: $GENOME"
        ((errors++))
    else
        echo "✓ Genome: $GENOME"
    fi
    
    if [ ! -f "$PEAKS" ]; then
        echo "✗ ERROR: Peaks file not found: $PEAKS"
        ((errors++))
    else
        echo "✓ Peaks: $PEAKS"
    fi
    
    if [ ! -f "$MOTIFS" ]; then
        echo "✗ ERROR: Motifs file not found: $MOTIFS"
        ((errors++))
    else
        echo "✓ Motifs: $MOTIFS"
    fi
    
    # Check optional files
    if [ -n "$BLACKLIST" ] && [ ! -f "$BLACKLIST" ]; then
        echo "⚠ WARNING: Blacklist file not found: $BLACKLIST"
        BLACKLIST=""
    elif [ -n "$BLACKLIST" ]; then
        echo "✓ Blacklist: $BLACKLIST"
    fi
    
    if [ -n "$PEAK_HEADER" ] && [ ! -f "$PEAK_HEADER" ]; then
        echo "⚠ WARNING: Peak header file not found: $PEAK_HEADER"
        PEAK_HEADER=""
    elif [ -n "$PEAK_HEADER" ]; then
        echo "✓ Peak header: $PEAK_HEADER"
    fi
    
    # Check samples
    if [ ${#SAMPLES[@]} -eq 0 ]; then
        echo "✗ ERROR: No samples specified"
        ((errors++))
    else
        echo "✓ Number of samples: ${#SAMPLES[@]}"
    fi
    
    # Check TOBIAS installation
    if ! command -v TOBIAS &> /dev/null; then
        echo "✗ ERROR: TOBIAS is not installed or not in PATH"
        ((errors++))
    else
        echo "✓ TOBIAS: $(which TOBIAS)"
    fi
    
    echo ""
    echo "Configuration Summary:"
    echo "  - Cores: $CORES"
    echo "  - Output directory: $OUTPUT_BASE"
    echo "  - Skip ATACorrect: $SKIP_ATACORRECT"
    echo "  - Skip FootprintScores: $SKIP_FOOTPRINT"
    echo "  - Skip BINDetect: $SKIP_BINDETECT"
    
    if [ $errors -gt 0 ]; then
        echo ""
        echo "✗ Validation failed with $errors error(s)"
        exit 1
    fi
    
    echo ""
    echo "✓ All validation checks passed!"
}

#############################################################
# Main Pipeline
#############################################################

run_pipeline() {
    print_header "Starting TOBIAS Pipeline"
    
    for sample_info in "${SAMPLES[@]}"; do
        # Parse sample information
        IFS='|' read -r SAMPLE_NAME BAM_FILE <<< "$sample_info"
        
        print_header "Processing Sample: $SAMPLE_NAME"
        
        # Check BAM file exists
        check_file "$BAM_FILE"
        
        # Define output directories
        OUTDIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
        SCORE_DIR="${OUTDIR}/score"
        BINDETECT_DIR="${OUTDIR}/BINDetect"
        
        # Create output directories
        mkdir -p "$OUTDIR"
        mkdir -p "$SCORE_DIR"
        mkdir -p "$BINDETECT_DIR"
        
        # Define output files
        CORRECTED_BW="${OUTDIR}/$(basename ${BAM_FILE%.bam})_corrected.bw"
        FOOTPRINT_BW="${SCORE_DIR}/${SAMPLE_NAME}_footprint.bw"
        
        #############################################################
        # Step 1: ATACorrect
        #############################################################
        
        if [ "$SKIP_ATACORRECT" = false ]; then
            print_step "Step 1/3: Running ATACorrect for $SAMPLE_NAME"
            
            if [ -f "$CORRECTED_BW" ]; then
                echo "⚠ Corrected bigWig already exists, skipping ATACorrect..."
                echo "  File: $CORRECTED_BW"
            else
                ATACORRECT_CMD="TOBIAS ATACorrect \
                    --bam \"$BAM_FILE\" \
                    --genome \"$GENOME\" \
                    --peaks \"$PEAKS\" \
                    --outdir \"$OUTDIR\" \
                    --cores $CORES"
                
                if [ -n "$BLACKLIST" ]; then
                    ATACORRECT_CMD="$ATACORRECT_CMD --blacklist \"$BLACKLIST\""
                fi
                
                eval $ATACORRECT_CMD
                
                echo "✓ ATACorrect completed for $SAMPLE_NAME"
            fi
        else
            echo "⊘ Skipping ATACorrect (--skip-atacorrect specified)"
            if [ ! -f "$CORRECTED_BW" ]; then
                echo "ERROR: Corrected bigWig not found: $CORRECTED_BW"
                exit 1
            fi
        fi
        
        #############################################################
        # Step 2: FootprintScores
        #############################################################
        
        if [ "$SKIP_FOOTPRINT" = false ]; then
            print_step "Step 2/3: Running FootprintScores for $SAMPLE_NAME"
            
            if [ -f "$FOOTPRINT_BW" ]; then
                echo "⚠ Footprint scores already exist, skipping FootprintScores..."
                echo "  File: $FOOTPRINT_BW"
            else
                TOBIAS FootprintScores \
                    --signal "$CORRECTED_BW" \
                    --regions "$PEAKS" \
                    --output "$FOOTPRINT_BW" \
                    --cores "$CORES"
                
                echo "✓ FootprintScores completed for $SAMPLE_NAME"
            fi
        else
            echo "⊘ Skipping FootprintScores (--skip-footprint specified)"
            if [ ! -f "$FOOTPRINT_BW" ]; then
                echo "ERROR: Footprint bigWig not found: $FOOTPRINT_BW"
                exit 1
            fi
        fi
        
        #############################################################
        # Step 3: BINDetect
        #############################################################
        
        if [ "$SKIP_BINDETECT" = false ]; then
            print_step "Step 3/3: Running BINDetect for $SAMPLE_NAME"
            
            if [ -d "$BINDETECT_DIR" ] && [ "$(ls -A $BINDETECT_DIR)" ]; then
                echo "⚠ BINDetect output directory exists and is not empty"
                echo "  Directory: $BINDETECT_DIR"
                echo "  Skipping BINDetect for $SAMPLE_NAME"
            else
                BINDETECT_CMD="TOBIAS BINDetect \
                    --motifs \"$MOTIFS\" \
                    --signals \"$FOOTPRINT_BW\" \
                    --genome \"$GENOME\" \
                    --peaks \"$PEAKS\" \
                    --outdir \"$BINDETECT_DIR\" \
                    --cond_names \"$SAMPLE_NAME\" \
                    --cores $CORES"
                
                if [ -n "$PEAK_HEADER" ]; then
                    BINDETECT_CMD="$BINDETECT_CMD --peak_header \"$PEAK_HEADER\""
                fi
                
                eval $BINDETECT_CMD
                
                echo "✓ BINDetect completed for $SAMPLE_NAME"
            fi
        else
            echo "⊘ Skipping BINDetect (--skip-bindetect specified)"
        fi
        
        print_step "Sample $SAMPLE_NAME completed successfully!"
        echo "Output directory: $OUTDIR"
        echo "  - Corrected bigWig: $CORRECTED_BW"
        echo "  - Footprint scores: $FOOTPRINT_BW"
        echo "  - BINDetect results: $BINDETECT_DIR"
        
    done
}

#############################################################
# Main Execution
#############################################################

main() {
    # Parse command line arguments first
    parse_arguments "$@"
    
    # Load configuration file if specified
    load_config
    
    # Load samples
    load_samples
    
    # Validate configuration
    validate_config
    
    # Run pipeline
    run_pipeline
    
    # Final summary
    print_header "TOBIAS Pipeline Completed Successfully!"
    echo "All samples processed: ${#SAMPLES[@]}"
    echo "Results saved in: $OUTPUT_BASE"
    echo ""
    echo "Next steps:"
    echo "  1. Check BINDetect results in each sample's BINDetect directory"
    echo "  2. Compare differential binding across samples"
    echo "  3. Visualize footprints for specific motifs"
    echo ""
}

# Run main function
main "$@"
