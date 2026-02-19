#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Footprint & SNP Integration Pipeline
================================================================================
Author:
    Thx
Date:
    2025-11-26
Description:
    This pipeline integrates ATAC-seq footprint scores with SNP (common and rare variants) annotations from
    dbSNP VCF files. It performs three main steps:
    1. BED Processing: 
       - Reads footprint scores in BED file from TOBIAS output(step3_filtered.bed from Tobias_filter_footprint_region.py)
       - Calculates relative scores using Min-Max normalization (0-1 scale)
       - Removes 'chr' prefix from chromosome names for VCFtools compatibility
       - Builds an in-memory lookup table for fast score retrieval
    2. VCF Intersection:
       - Uses VCFtools to extract SNPs that overlap with footprint regions
       - Filters out indels, keeping only SNPs
       - Preserves all INFO field annotations
    3. Annotation Merging:
       - Parses filtered VCF file
       - Extracts comprehensive SNP annotations (gene info, frequencies, etc.)
       - Appends relative footprint scores to each SNP
       - Outputs final CSV with integrated data
Usage:
    Basic usage:
        python footprint_snp_pipeline.py \\
            --bed footprints.bed \\
            --vcf dbsnp151.vcf \\
            --out-prefix output_nam
    Example:
        python footprint_snp_pipeline.py \\
            --bed HEPM_ATAC_1_Q75_step3_filtered.bed \\
            --vcf dbsnp151.vcf \\
            --out-prefix AEs_805_TFbSNP
Parameters:
    Required Arguments:
        --bed FILE          Input BED file with footprint scores
                           Format: chr start end score [optional_columns]
                           Example: chr1 3192243 3192244 57.45974
        
        --vcf FILE          Input dbSNP VCF file (can be bgzipped)
                           Must contain INFO fields: GENEINFO, TOPMED, etc.
                           Example: dbsnp151.vcf or dbsnp151.vcf.gz
    Optional Arguments:
        --out-prefix STR    Output file prefix (default: AEs_805_TFbSNP)
                           Generates multiple output files:
                           - {prefix}.csv: Final annotated SNPs
                           - {prefix}.recode.vcf: Filtered VCF
                           - {prefix}_nochr.bed: Processed BED (temp)
Output Files:
    1. {prefix}.csv
       - Main output: Annotated SNPs with footprint scores
       - Columns: CHROM, POS, ID, REF, ALT, GENE_SYMBOL, GENE_ID,
                 VARIANT_TYPE, REF_FREQ, ALT_FREQ, FREQUENCY_CLASS,
                 RELATIVE_FOOTPRINT_SCORE, HD, GNO, KGPhase1, KGPhase3, SSR
    2. {prefix}.recode.vcf
    3. {bed_name}_nochr.bed
Notes:
    - BED coordinates are 0-based, half-open [start, end)
    - VCF coordinates are 1-based, inclusive
    - SNP matching: BED region [start, end) contains VCF POS if start < POS <= end
    - Relative scores are Min-Max normalized to [0, 1] range
    - SNPs with SSR > 0 (suspect regions) are filtered out

"""

import csv
import argparse
import subprocess
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
from collections import Counter, defaultdict

# ==========================================
# Part 1: BED Processing Class
# ==========================================

class BEDProcessor:
    """Handles BED file normalization and formatting"""
    
    def __init__(self, input_bed: str):
        self.input_bed = input_bed
        self.score_map = defaultdict(list)
        
    def process_and_save(self, output_bed: str) -> None:
        """
        Reads input BED, calculates relative score (Min-Max scaling),
        removes 'chr' prefix, and saves to output_bed.
        Also populates self.score_map for later lookup.
        """
        print(f"🔄 Processing BED file: {self.input_bed}...")
        
        raw_data = []
        scores = []
        
        # 1. Read Data
        try:
            with open(self.input_bed, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) < 4:
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    score = float(parts[3])
                    
                    raw_data.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'score': score,
                        'rest': parts[4:] # Keep other columns if any
                    })
                    scores.append(score)
        except FileNotFoundError:
            print(f"❌ Error: BED file not found: {self.input_bed}")
            sys.exit(1)

        if not scores:
            print("❌ Error: No valid data found in BED file.")
            sys.exit(1)

        # 2. Calculate Relative Score (Min-Max Normalization: 0 to 1)
        min_s = min(scores)
        max_s = max(scores)
        range_s = max_s - min_s if max_s != min_s else 1.0
        
        print(f"   📊 Score Range: Min={min_s}, Max={max_s}")
        
        # 3. Write processed BED and build Lookup Dict
        with open(output_bed, 'w') as f_out:
            for item in raw_data:
                # Calculate relative score
                rel_score = (item['score'] - min_s) / range_s
                
                # Remove 'chr' prefix for VCFtools compatibility
                chrom_nochr = item['chrom'].replace('chr', '')
                
                # Store in memory for VCF merging later
                # Structure: { '1': [(3192243, 3192244, 0.55), ...], '2': ... }
                self.score_map[chrom_nochr].append((item['start'], item['end'], rel_score))
                
                # Write to file: Chrom Start End Score RelativeScore
                # Note: We write the relative score as the 5th column
                line_out = f"{chrom_nochr}\t{item['start']}\t{item['end']}\t{item['score']}\t{rel_score:.6f}\n"
                f_out.write(line_out)
                
        print(f"✅ Processed BED saved to: {output_bed}")

    def get_score_for_position(self, chrom_nochr: str, pos: int) -> float:
        """
        Finds the relative score for a specific SNP position.
        Checks if Start < POS <= End (Standard BED to VCF mapping).
        """
        # Optimization: If the list is sorted, we could use binary search.
        # For simplicity and safety with overlapping regions, we iterate.
        # Given VCFtools already filtered, the SNP *must* be in one of these regions.
        
        if chrom_nochr not in self.score_map:
            return 0.0
            
        for start, end, rel_score in self.score_map[chrom_nochr]:
            # VCF POS is 1-based. BED Start is 0-based, End is 1-based.
            # A SNP at 100 falls in BED [99, 100).
            # Usually for footprinting: BED End == VCF POS.
            # Let's use standard interval check:
            if start < pos <= end:
                return rel_score
        
        return 0.0 # Should not happen if VCFtools worked correctly

# ==========================================
# Part 2: VCF Parsing Class (Modified)
# ==========================================

@dataclass
class SNPRecord:
    """Data class to store SNP information"""
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    GENE_SYMBOL: str
    GENE_ID: str
    VARIANT_TYPE: str
    REF_FREQ: float
    ALT_FREQ: float
    FREQUENCY_CLASS: str
    HD: bool
    GNO: bool
    KGPhase1: bool
    KGPhase3: bool
    SSR: int
    RELATIVE_FOOTPRINT_SCORE: float  # <--- Added this field
    
    def to_dict(self) -> Dict:
        return asdict(self)

class VCFParser:
    """Parser for VCF files with SNP annotation extraction"""
    
    def __init__(self, bed_processor: BEDProcessor, rare_threshold: float = 0.01):
        self.rare_threshold = rare_threshold
        self.bed_processor = bed_processor # Link to BED data
        self.stats = {
            'total_variants': 0,
            'filtered_by_ssr': 0,
            'variant_types': Counter(),
            'frequency_classes': Counter(),
            'missing_topmed': 0,
            'missing_geneinfo': 0
        }
    
    # ... [Previous helper methods: parse_info_field, extract_geneinfo, etc. remain same] ...
    def parse_info_field(self, info_str: str) -> Dict[str, str]:
        info_dict = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
            else:
                info_dict[item] = 'True'
        return info_dict
    
    def extract_geneinfo(self, info_dict: Dict[str, str]) -> Tuple[str, str]:
        geneinfo = info_dict.get('GENEINFO', '')
        if not geneinfo:
            self.stats['missing_geneinfo'] += 1
            return ('', '')
        genes = geneinfo.split('|')
        if ':' in genes[0]:
            gene_symbol, gene_id = genes[0].split(':', 1)
            return (gene_symbol, gene_id)
        return (geneinfo, '')
    
    def extract_variant_type(self, info_dict: Dict[str, str]) -> str:
        variant_types = []
        for t in ['INT', 'SYN', 'NSN', 'NSM', 'U3', 'U5', 'R3', 'R5']:
            if t in info_dict: variant_types.append(t)
        if variant_types:
            variant_type = '|'.join(variant_types)
        else:
            variant_type = info_dict.get('VC', 'UNKNOWN')
        self.stats['variant_types'][variant_type] += 1
        return variant_type
    
    def extract_topmed_freq(self, info_dict: Dict[str, str]) -> Tuple[float, float]:
        topmed = info_dict.get('TOPMED', '')
        if not topmed:
            self.stats['missing_topmed'] += 1
            return (0.0, 0.0)
        freqs = topmed.split(',')
        try:
            ref_freq = float(freqs[0])
            alt_freq = float(freqs[1]) if len(freqs) > 1 else 0.0
            return (ref_freq, alt_freq)
        except (ValueError, IndexError):
            self.stats['missing_topmed'] += 1
            return (0.0, 0.0)
    
    def classify_frequency(self, alt_freq: float) -> str:
        if alt_freq == 0.0: freq_class = 'UNKNOWN'
        elif alt_freq < self.rare_threshold: freq_class = 'RARE'
        else: freq_class = 'COMMON'
        self.stats['frequency_classes'][freq_class] += 1
        return freq_class
    
    def check_boolean_flags(self, info_dict: Dict[str, str]) -> Dict[str, bool]:
        return {
            'HD': 'HD' in info_dict, 'GNO': 'GNO' in info_dict,
            'KGPhase1': 'KGPhase1' in info_dict, 'KGPhase3': 'KGPhase3' in info_dict
        }

    def parse_vcf_line(self, line: str) -> Optional[SNPRecord]:
        fields = line.strip().split('\t')
        if len(fields) < 8: return None
        
        chrom, pos, snp_id, ref, alt = fields[0:5]
        info_str = fields[7]
        
        info_dict = self.parse_info_field(info_str)
        
        ssr = int(info_dict.get('SSR', '0'))
        if ssr > 0:
            self.stats['filtered_by_ssr'] += 1
            return None
        
        # Handle Chromosome name for output (add 'chr' back for final CSV if desired)
        # But for lookup in BED map, we need it without 'chr' (as we stripped it in BEDProcessor)
        chrom_lookup = chrom.replace('chr', '')
        pos_int = int(pos)
        
        # *** Retrieve Relative Score from BED Processor ***
        rel_score = self.bed_processor.get_score_for_position(chrom_lookup, pos_int)
        
        # Add 'chr' prefix for final output consistency
        if not chrom.startswith('chr'):
            chrom_out = f'chr{chrom}'
        else:
            chrom_out = chrom

        gene_symbol, gene_id = self.extract_geneinfo(info_dict)
        variant_type = self.extract_variant_type(info_dict)
        ref_freq, alt_freq = self.extract_topmed_freq(info_dict)
        freq_class = self.classify_frequency(alt_freq)
        boolean_flags = self.check_boolean_flags(info_dict)
        
        record = SNPRecord(
            CHROM=chrom_out,
            POS=pos_int,
            ID=snp_id,
            REF=ref,
            ALT=alt,
            GENE_SYMBOL=gene_symbol,
            GENE_ID=gene_id,
            VARIANT_TYPE=variant_type,
            REF_FREQ=ref_freq,
            ALT_FREQ=alt_freq,
            FREQUENCY_CLASS=freq_class,
            HD=boolean_flags['HD'],
            GNO=boolean_flags['GNO'],
            KGPhase1=boolean_flags['KGPhase1'],
            KGPhase3=boolean_flags['KGPhase3'],
            SSR=ssr,
            RELATIVE_FOOTPRINT_SCORE=rel_score # <--- Set the score
        )
        
        self.stats['total_variants'] += 1
        return record
    
    def parse_vcf_file(self, vcf_file: str) -> List[SNPRecord]:
        records = []
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                record = self.parse_vcf_line(line)
                if record: records.append(record)
        return records
    
    def write_csv(self, records: List[SNPRecord], output_file: str):
        if not records:
            print("Warning: No records to write!")
            return
        
        fieldnames = [
            'CHROM', 'POS', 'ID', 'REF', 'ALT',
            'GENE_SYMBOL', 'GENE_ID', 'VARIANT_TYPE',
            'REF_FREQ', 'ALT_FREQ', 'FREQUENCY_CLASS',
            'RELATIVE_FOOTPRINT_SCORE', 
            'HD', 'GNO', 'KGPhase1', 'KGPhase3', 'SSR'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for record in records:
                writer.writerow(record.to_dict())
        print(f"✅ Successfully wrote {len(records)} records to {output_file}")

# ==========================================
# Part 3: Main Execution Pipeline
# ==========================================

def main():
    parser = argparse.ArgumentParser(description='Integrated Footprint & SNP Analysis Pipeline')
    
    # Inputs
    parser.add_argument('--bed', required=True, help='Input BED file (with footprint scores)')
    parser.add_argument('--vcf', required=True, help='Input dbSNP VCF file')
    
    # Outputs
    parser.add_argument('--out-prefix', default='AEs_805_TFbSNP', help='Output prefix for files')
    
    args = parser.parse_args()
    
    # File Names
    temp_bed_nochr = f"{args.bed.replace('.bed', '')}_nochr.bed"
    vcf_out_prefix = args.out_prefix # vcftools adds .recode.vcf
    final_csv = f"{args.out_prefix}.csv"
    
    # --- Step 1: Process BED ---
    print("\n🔹 STEP 1: Processing BED File...")
    bed_proc = BEDProcessor(args.bed)
    bed_proc.process_and_save(temp_bed_nochr)
    
    # --- Step 2: Run VCFtools ---
    print("\n🔹 STEP 2: Running VCFtools...")
    # Note: vcftools automatically appends ".recode.vcf" to the --out parameter
    vcftools_cmd = [
        "vcftools",
        "--vcf", args.vcf,
        "--bed", temp_bed_nochr,
        "--remove-indels",
        "--recode",
        "--recode-INFO-all",
        "--out", vcf_out_prefix
    ]
    
    cmd_str = " ".join(vcftools_cmd)
    print(f"   Command: {cmd_str}")
    
    try:
        subprocess.run(vcftools_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ VCFtools failed with error: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("❌ Error: vcftools not found in PATH.")
        sys.exit(1)
        
    intersected_vcf = f"{vcf_out_prefix}.recode.vcf"
    if not os.path.exists(intersected_vcf):
        print(f"❌ Error: Expected output VCF {intersected_vcf} not found.")
        sys.exit(1)
        
    print(f"✅ VCF Intersection complete. Output: {intersected_vcf}")
    
    # --- Step 3: Parse VCF & Merge Scores ---
    print("\n🔹 STEP 3: Parsing VCF and Merging Scores...")
    vcf_parser = VCFParser(bed_processor=bed_proc) # Pass the bed processor to access scores
    
    records = vcf_parser.parse_vcf_file(intersected_vcf)
    vcf_parser.write_csv(records, final_csv)
    
    # Cleanup (Optional)
    # os.remove(temp_bed_nochr)
    # os.remove(intersected_vcf)
    
    print("\n🎉 Pipeline Completed Successfully!")

if __name__ == '__main__':
    main()
