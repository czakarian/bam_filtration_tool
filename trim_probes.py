#!/usr/bin/env python

import re 
import pysam
import Bioinfo

class ProbeTrim():
    """The ProbeTrim class is used to trim off probe sequences (and corresponding quality scores) from alignments in a bam file. 
    The trimmed probe sequence will be appended to the header. Cigar string and start positions will also be adjusted accordingly."""

    def __init__(self, line, probe_length=40, side=3):
        """Initialize a ProbeTrim object with the following:
        line - pysam AlignmentSegment object, 
        probe_length - expected probe sequence length (default=40),
        side - whether probe is found at 5' or 3' end of sequence (default=3)."""
        self.line = line
        self.probe_length = probe_length
        self.side = side 

        # determine amount to trim
        self.trim = self.get_trim_length()
 
        # if the actual sequence is shorter than the amount to trim return None
        if self.trim >= len(self.line.query_sequence):
            self.line = None
        else:
            self.trim_seq_and_qs()
            self.trim_start()
            self.trim_cigar()
            self.trim_template_length()

    # def is_valid_alignment(self):
    #     if len(self.line.query_sequence) == 0:
    #         return None

    def get_trim_length(self):
        """This function returns the adjusted sequence length to trim off if there are any insertions
        or deletions encountered within the expected probe region."""

        trim_amount = self.probe_length

        # get a list of each of the different cigar string components
        result = re.findall('([0-9]+)([MDNIS])', self.line.cigarstring)

        # reverse cigar for 3' for or 5' rev strands
        if (self.side==3 and not self.line.is_reverse) or (self.side==5 and self.line.is_reverse):
            result.reverse()

        # determine how much probe to trim (if any insertions or deletions 40 +/- ?)
        sum = 0    
        for match in result:
            if sum + int(match[0]) <= self.probe_length:
                if match[1] == 'I':
                   trim_amount += int(match[0])
                if match[1] == 'D':
                   trim_amount -= int(match[0])
            sum += int(match[0]) 
        return trim_amount

    def trim_seq_and_qs(self):
        """This function trims the sequence and quality score lines and appends the trimmed sequence to the header."""
        seq = self.line.query_sequence
        qs = self.line.query_qualities
        # trim the seq and qs lines + append trimmed seq to header
        if (self.side==3 and not self.line.is_reverse) or (self.side==5 and self.line.is_reverse):
            self.line.query_sequence = seq[0:-self.trim]
            self.line.query_qualities = qs[0:-self.trim]
            self.line.query_name = self.line.query_name + '_' + Bioinfo.rev_comp(seq[-self.trim:])
        elif (self.side==3 and self.line.is_reverse) or (self.side==5 and not self.line.is_reverse):
            self.line.query_sequence = seq[self.trim:]
            self.line.query_qualities = qs[self.trim:]
            self.line.query_name = self.line.query_name + '_' + seq[0:self.trim]

    def trim_start(self):
        """This function adjusts the alignment start position based on amount of sequence to be trimmed."""
        if (self.side==3 and self.line.is_reverse) or (self.side==5 and not self.line.is_reverse):
            self.line.reference_start = self.line.reference_start + self.probe_length

    def trim_template_length(self):
        """This function adjusts the template length (for PE reads) based on amount of sequence to be trimmed."""
        if self.line.template_length < 0:
            self.line.template_length = self.line.template_length + self.trim
        elif self.line.template_length > 0:
            self.line.template_length = self.line.template_length - self.trim
        else:
            # if tem_len is 0, not paired end
            pass

    def trim_cigar(self) :
        """This function adjusts the cigar string based on amount of sequence to be trimmed. """

        temp = 0
        first = True
        new_cigar = []
        # get a list of each of the different cigar string components
        result = re.findall('([0-9]+)([MDNIS])', self.line.cigarstring)
        # temp reverse cigar for 3' for or 5' rev strands
        if (self.side==3 and not self.line.is_reverse) or (self.side==5 and self.line.is_reverse):
            result.reverse()

        for match in result:
            if temp + int(match[0]) > self.trim:
                if first:
                    if match[1] in {'I', 'N'}:  # if the new first cigar componenet is I or N, change it to soft-clipping
                        match = (match[0], 'S')
                    if match[1] == 'D':     # if new first cigar componenent is D, don't include it in the new cigar string
                        continue 
                    new_cigar.append(str(int(match[0]) - self.trim + temp) + match[1])
                    first = False
                else:
                    new_cigar.append(match[0] + match[1])
            # don't count deletions that have been removed as part of temp sum (since no actual bases were removed from seq)
            if match[1] != 'D':
                temp += int(match[0])

        # re-reverse cigar for 3' forward or 5' reverse stands
        if (self.side==3 and not self.line.is_reverse) or (self.side==5 and self.line.is_reverse):
            new_cigar.reverse()

        self.line.cigarstring = ''.join(new_cigar)
