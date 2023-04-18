This provided set of scripts is used to perform filtering of alignments in a bam file as well as alignment-based probe trimming. Alignments will be filtered based on whether they are 1. unmapped/not in proper pair, 2. supplementary 3. short (<65N), 4. off-target, 5. on-target. A separate bam file will be produced as output for each of these 5 categories.

Input 
argparse options:
    -i, --input: required arg, file path of BAM file
    -o, --out: required arg, file path to output directory (must exist)
    -t, --probes: required arg, file path to tsv file containing probe sequences used (column 1 = names, column 2 = sequences)
    -b, --bed: required arg, file path to bed file generated from probe sequences
    -p, --paired: optional arg, If flag is set, will assume given a paired end bam.

Output
5 BAM files as follows:
    {file}.on.bam // on-target 
    {file}.off.bam // off-target 
    {file}.sup.bam // supplementary 
    {file}.len.bam // short length (<65N)
    {file}.unmap.bam // unmapped or not in proper pair
2 stats files:
    {file}.lev_dists.tsv // Distribution of detected mismatches in probe sequences according to levenshtein distance
    {file}.trimmed_probe_lengths.tsv // Distribution of trimmed probe lengths (will signify presense of ins/del in probe seqs)

Will also output count data (# reads in each filtration category) to stdout.

Example command execution:
```
python run_filter_alignments.py -i sample_1.bam -o output_dir/ -t probes.tsv -b probes.bed --paired
```