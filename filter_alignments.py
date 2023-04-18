#!/usr/bin/env python

import re
import pysam
from Levenshtein import distance as lev
from trim_probes import ProbeTrim

class FilterAlignments():
    """The FilterAlignments class is used to filter alignments in a bam file (unmapped or not in proper pair, supplementary, short alignments(<65N), off-target, and on-target reads. 
    It will output seperate bam files containing the alignments that belong to each of the mentioned categories."""

    def __init__(self, **kwargs):
        self.raw_input_seq_path = kwargs.get('input_seq_path')
        self.raw_output_dir = kwargs.get('output_dir')
        self.raw_probes_tsv_path = kwargs.get('probe_tsv_path')
        self.raw_probes_bed_path = kwargs.get('probe_bed_path')
        self.raw_paired= kwargs.get('paired_bam')

        # initialize dictionaries and counters
        self.lev_dists = {}
        self.trimmed_lengths = {} 
        self.sup_headers = set()
        self.pp_headers = set()
        self.counts = {'total':0, 'on_target':0, 'off_target':0, 'short':0, 'unmap':0, 'suppl':0, 'suppl_on':0, 'suppl_off':0}
        self.pair_line = ''

        # set cutoffs
        self.max_lev_dist = 3
        self.min_match_seq = 10
        self.min_seq_len = 65

        # read in probe tsv and probe bed 
        self.probes = self.read_probe_tsv(self.raw_probes_tsv_path)
        self.probe_positions = self.read_probe_bed(self.raw_probes_bed_path)

        self.file_prefix = '.'.join(self.raw_input_seq_path.split('/')[-1].split('.')[:-1])

        # open input and output bam files
        self.fr = pysam.AlignmentFile(self.raw_input_seq_path, 'rb')
        self.fw_on = pysam.AlignmentFile(self.raw_output_dir  + self.file_prefix + '.on.bam', 'wb', template = self.fr)
        self.fw_off = pysam.AlignmentFile(self.raw_output_dir  + self.file_prefix + '.off.bam', 'wb', template = self.fr)
        self.fw_sup = pysam.AlignmentFile(self.raw_output_dir  + self.file_prefix + '.sup.bam', 'wb', template = self.fr)
        self.fw_len = pysam.AlignmentFile(self.raw_output_dir  + self.file_prefix + '.len.bam', 'wb', template = self.fr)
        self.fw_unmap = pysam.AlignmentFile(self.raw_output_dir  + self.file_prefix + '.unmap.bam', 'wb', template = self.fr)

        self.parse_file()
        self.write_stats_to_file()

        self.fr.close()
        self.fw_on.close() 
        self.fw_off.close() 
        self.fw_sup.close()
        self.fw_len.close()
        self.fw_unmap.close()

    def get_5prime(self, start_pos:int, cigar:str, reverse_strand:bool) -> int:
        """This function returns the true 5' start position (adjusted for soft clipping) for reads mapping to + strand or the - strand."""
        # for + strand reads, adjust soft clipping at the beginning of CIGAR string
        if not reverse_strand:
            result = re.match('^([0-9]+)S', cigar)
            if result:  
                start_pos -= int(result[1])  
        # for - strand reads, sum M/D/N/S values and add that to start_pos, ignore any inserts and soft clipping at the beginning   
        else:
            result_S = re.search('([0-9]+)S$', cigar)
            result_MDN = re.findall('([0-9]+)[MDN]', cigar)
            if result_S:
                start_pos += int(result_S[1])
            if result_MDN:
                start_pos += sum(map(int, result_MDN))
            # subtract 1 to correct to actual position 
        start_pos -= 1  
        return start_pos

    def get_3prime(self, start_pos:int, cigar:str, reverse_strand:bool) -> int:
        """This function returns the 3' end position (adjusted for soft clipping) for reads mapping to + strand or the - strand."""
        # get end_pos from start_pos
        end_pos = start_pos
        # for - strand reads, adjust soft clipping at the beginning of CIGAR string
        if reverse_strand:
            result = re.match('^([0-9]+)S', cigar)
            if result:  
                end_pos -= int(result[1])
        # for + strand reads, sum M/D/N/S values and add that to start_pos, ignore any inserts and soft clipping at the beginning   
        else:
            result_S = re.search('([0-9]+)S$', cigar)
            result_MDN = re.findall('([0-9]+)[MDN]', cigar)
            if result_S:
                end_pos += int(result_S[1])
            if result_MDN:
                end_pos += sum(map(int, result_MDN))
            # subtract 1 to correct to actual position
        end_pos -= 1
        return end_pos

    def pop_freq_dict(self,item, my_dict):
        if item in my_dict:
            my_dict[item] += 1
        else:
            my_dict[item] = 1

    def get_probe_lev_dist(self, line):
        # get levenstein distance between trimmed probe and error-corrected probe
        trimmed_probe = line.query_name.split('_')[-2]
        ec_probe = line.query_name.split('_')[-1]
        ld = lev(trimmed_probe, ec_probe)

        # add lev dist and trimmed length to count dictionaries
        self.pop_freq_dict(ld, self.lev_dists)
        self.pop_freq_dict(len(trimmed_probe), self.trimmed_lengths)

        return ld

    def read_probe_bed(self, file):
        """read in probe bed file
        file structure
            1	900053	900093	chr1_rs4422948_f	60	+
            1	900144	900184	chr1_rs4422948_r	60	-
        """
        probe_dict = {}
        with open(file, 'r') as fr:  
            for line in fr:
                cols = line.strip().split()
                if cols[5] == '+':
                    probe_dict[(cols[0], cols[1], cols[5])] = cols[3]
                elif cols[5] == '-':
                    probe_dict[(cols[0], cols[2], cols[5])] = cols[3]
        return probe_dict

    def read_probe_tsv(self, file):
        """read in probe tsv file
        file structure	
            ALKintron19_1	CTGCAGCTCCATCTGCATGGCTTGCAGCTCCTGGTGCTTC
            ALKintron19_2	CGGCGGTACACTGCAGGTGGGTGGTCAGCTGCAACATGGC
        """
        probe_dict = {}
        with open(file, 'r') as fr:   
            fr.readline()   
            for line in fr:  
                cols = line.strip().split()
                probe_dict[cols[0]] = cols[1].upper()
            return probe_dict

    def is_properly_mapping_align(self, line):
        # Filter out reads that are not in proper pair
        if not line.is_proper_pair:
            self.fw_unmap.write(line)
            if line.query_name not in self.pp_headers:
                self.counts['unmap'] += 1
                self.pp_headers.add(line.query_name)                
            return False 
        return True

    def is_suppl_align(self, line):
        # Filter out supplementary alignments
        if line.has_tag('SA'):
            if line.is_read2 and self.pair_line:
                self.fw_sup.write(self.pair_line)
                self.pair_line = None
            self.fw_sup.write(line)
            if line.query_name not in self.sup_headers:
                self.counts['suppl'] += 1
                self.sup_headers.add(line.query_name)
            return True
        # if alignment not marked with SA:Z, but header has been seen, alignment is a read 2 pair to a read 1 suppl. alignment
        elif not line.has_tag('SA') and line.query_name in self.sup_headers:
            self.fw_sup.write(line)
            return True    
        return False      

    def is_short_align(self, line):
        # Filter any alignments <65N (assuming ~40N should be probe sequence)
        if abs(line.tlen) <= 65:
            #self.fw_len.write(self.pair_line)
            #self.pair_line = None
            self.fw_len.write(line)
            if line.is_read2:
                self.counts['short'] += 1
            return True
        return False
    
    def parse_file(self):
        for line in self.fr: 
            probe_strand = ''

            if self.raw_paired:
                if not self.is_properly_mapping_align(line) or self.is_suppl_align(line) or self.is_short_align(line):
                    continue

                # if read 1 is not part of a suppl, temp store until next iteration (read 2) to determine whether read 2 is suppl/on/off and write stored read 1 to appropriate file 
                if line.is_read1:
                    self.pair_line = line
                    continue 

                # Determine 5' position of read 2 to get probe end position
                probe_end = self.get_5prime(line.reference_start + 1, line.cigarstring, line.is_reverse)

                # if read 2 is reverse mapping, then the actual probe sequence should map to reverse strand and vice versa
                if line.is_reverse:
                    probe_strand = '-'
                else:
                    probe_strand = '+'

                # If 5' position of read 2 maps to a known probe end position, proceed
                result = None
                if (str(line.reference_name), str(probe_end), probe_strand) in self.probe_positions:
                    # probe trim will trim probe seq and append to the trimmed sequence to header
                    trim_line = ProbeTrim(line, side=5).line

                    # Filter any alignments that get completely trimmed
                    if trim_line == None:
                        self.counts['short'] += 1
                        continue   

                    # also append correct probe sequence to header (ex. header_tag_trimprobe_correctprobe)
                    trim_line.query_name = trim_line.query_name + '_' + self.probes[self.probe_positions[(str(trim_line.reference_name), str(probe_end), probe_strand)]]
                    # add in the trimmed and ec probe sequences to read 1 of the read pair as well
                    self.pair_line.query_name = trim_line.query_name

                    # adjust template length for read 1, according to adjusted template length of read 2 
                    self.pair_line.template_length = -trim_line.template_length
                    
                    if trim_line.is_reverse:
                        result = re.search('([0-9]+)M$', trim_line.cigarstring)
                    else:
                        result = re.match('^([0-9]+)M', trim_line.cigarstring)

                # require at least a 10M sequence before the trimmed probe, otherwise mark as off target
                if result and int(result[1]) >= self.min_match_seq:

                    # get lev dist between trimmed probe and ec probe
                    ld = self.get_probe_lev_dist(trim_line)

                    # allow a max lev distance of 3, otherwise mark as off target
                    if ld <= self.max_lev_dist:
                        self.counts['on_target'] += 1
                        self.fw_on.write(self.pair_line)
                        self.fw_on.write(trim_line)
                    else:
                        self.counts['off_target'] += 1
                        self.fw_off.write(self.pair_line)
                        self.fw_off.write(line)

                else:
                    self.counts['off_target'] += 1
                    self.fw_off.write(self.pair_line)
                    self.fw_off.write(line)

            else:
                # Filter out unmapped 
                if line.is_unmapped:
                    self.counts['unmap'] += 1
                    self.fw_unmap.write(line)
                    continue

                # Filter out supplementary / chimeric alignments        
                if line.has_tag('SA'):
                    if line.query_name not in self.sup_headers:
                        self.counts['suppl'] += 1
                        self.sup_headers.add(line.query_name)
                    self.fw_sup.write(line)
                    continue

                # Filter any alignments <65N (assuming ~40N should be probe sequence)
                if len(line.query_sequence) <= 65:
                    self.fw_len.write(line)
                    self.counts['short'] += 1
                    continue

                # Determine 3' end position to get probe end position
                probe_end = self.get_3prime(line.reference_start + 1, line.cigarstring, line.is_reverse)

                # if merged read is forward mapping, then the actual probe sequence should map to reverse strand and vice versa
                if not line.is_reverse:
                    probe_strand = '-'
                else:
                    probe_strand = '+'

                # if 3' end of read maps to a known probe end position, proceed
                result = None
                if (str(line.reference_name), str(probe_end), probe_strand) in self.probe_positions:
                    # probe trim will trim probe seq and append trimmed sequence to header
                    trim_line = ProbeTrim(line, side=3).line

                    # also append correct probe sequence to header (ex. header_tag_trimprobe_correctprobe)
                    trim_line.query_name = trim_line.query_name + '_' + self.probes[self.probe_positions[(str(trim_line.reference_name), str(probe_end), probe_strand)]]
                    
                    if not trim_line.is_reverse:
                        result = re.search('([0-9]+)M$', trim_line.cigarstring)
                    else:
                        result = re.match('^([0-9]+)M', trim_line.cigarstring)

                # require at least a 10M sequence before the trimmed probe, otherwise mark as off target
                if result and int(result[1]) >= self.min_match_seq:

                    # get lev dist between trimmed probe and ec probe
                    ld = self.get_probe_lev_dist(trim_line)

                    # allow a max lev distance of 3, otherwise mark as off target
                    if ld <= self.max_lev_dist:
                        self.counts['on_target'] += 1
                        self.fw_on.write(trim_line)
                    else:
                        self.counts['off_target'] += 1
                        self.fw_off.write(line)
                else:
                    self.counts['off_target'] += 1
                    self.fw_off.write(line)

    def write_stats_to_file(self):
        with open(self.raw_output_dir + self.file_prefix + '.lev_dists.tsv', 'w') as fw:
            fw.write('Distance' + '\t' + 'Frequency' + '\n') 
            for i in sorted(self.lev_dists):
                fw.write(str(i) + '\t' + str(self.lev_dists[i]) + '\n')

        with open(self.raw_output_dir + self.file_prefix + '.trimmed_probe_lengths.tsv', 'w') as fw:
            fw.write('Length' + '\t' + 'Frequency' + '\n') 
            for i in sorted(self.trimmed_lengths):
                fw.write(str(i) + '\t' + str(self.trimmed_lengths[i]) + '\n')

        self.counts['total'] = self.counts['off_target'] + self.counts['on_target'] + self.counts['suppl']  + self.counts['short']  + self.counts['unmap'] 
        print('Total Reads: ', self.counts['total'])
        print('Unmapped Reads (or not in proper pair for PE): ', self.counts['unmap'], ',', round(self.counts['unmap'] / self.counts['total'] * 100,2))
        print('Supplementary Reads: ', self.counts['suppl'], ',', round(self.counts['suppl'] / self.counts['total'] * 100,2))
        print('Short Reads (<65N): ', self.counts['short'], ',', round(self.counts['short'] / self.counts['total'] * 100,2))
        if(self.counts['on_target'] or self.counts['off_target']):
            print('On-target Reads: ', self.counts['on_target'], ',', round((self.counts['on_target'] / (self.counts['on_target'] + self.counts['off_target'])) * 100,2))
            print('Off-target Reads: ', self.counts['off_target'], ',', round((self.counts['off_target'] / (self.counts['on_target'] + self.counts['off_target'])) * 100,2))

