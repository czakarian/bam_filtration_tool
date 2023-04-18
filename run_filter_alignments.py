#!/usr/bin/env python

"""This script takes in a bam file (specify with flag if paired end) and filters alignments that are off-target, supplementary, 
unmapped / not in proper pair, or short alignments (<65N). Remaining on-target reads will have probes trimmed before being written to output file.
An alignment will be considered on-target according to alignment positions in the provided probe bed file, if it has at least 10N matching seq
before the probe, and has 3 or fewer mismatches in the probe sequence. The input bam should not have had probes trimmed as this script assumes
the probe is present in the alignment."""

import os
import sys
import argparse

from filter_alignments import FilterAlignments

def get_args(argv=None):
    parser = argparse.ArgumentParser(description="This program performs filtering of off-target alignments, supplementary alignments, unmapped / not in proper pair, length filtering, and probe trimming.")
    parser.add_argument("-i", "--input", help="Input bam file.", required=True)
    parser.add_argument("-o", "--out", help="Path to output directory.", required=True)
    parser.add_argument("-t", "--probes", help="Path to probe tsv.", required=True)
    parser.add_argument("-b", "--bed", help="Path to probe bed file.", required=True)
    parser.add_argument("-p", "--paired", help="If flag is set, will assume given a paired end bam.", action='store_true')    
    return parser.parse_args()

def filter_alignments_main(argv=None):
    args = get_args(argv)

    if args.out[-1] != '/':
        args.out = args.out + '/'

    filter_args = {'input_seq_path': args.input,
                  'output_dir': args.out,
                  'probe_tsv_path': args.probes,
                  'probe_bed_path': args.bed,
                  'paired_bam': args.paired}

    filtration_obj = FilterAlignments(**filter_args)
    #filtration_obj.run_analysis()

if __name__ == "__main__":
    filter_alignments_main()