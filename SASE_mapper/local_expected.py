import numpy as np
from pybedtools import BedTool
from collections import defaultdict
import time


def calc_expected(evr_lengths, evr_sums, long_overlap):
    '''
    Given the EVR proportions, calculate expected number of
    mutations in potential SASE
    
    requires
        evr_lengths: dictionary od EVR labels -> total length of EVR
        evr_sums: dictionary of EVR labels -> sum of p-values occurring in EVR
        long_overlap: list of dictionaries for number of basepairs overlapping EVR
    return
        local_expected: array of expected values
    '''
    
    evr_expected = {}    
    for label in evr_sums:
        evr_expected[label] = evr_sums[label] / float(evr_lengths[label])
    
    #print evr_expected
    
    local_expected = np.zeros(len(long_overlap), dtype=float)
    for i in xrange(len(long_overlap)):
        expected = 0

        for label in long_overlap[i]:
            expected += evr_expected[label] * long_overlap[i][label]
    
        local_expected[i] = expected

    return local_expected


def generate_long(long_fn):
    '''
    requires
        long_fn: long output file name (contains potential SASEs)
    yields: tuple(chrom name, start, end, signal)
    '''
    
    for line in open(long_fn, "r"):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        pvalue = float(fields[3])
        
        if pvalue < 1.0:
            signal = -np.log10(max(1e-50, pvalue))
        
            yield chrom, start, end, signal


def calculate_local_expected(long_fn, EVR_fn):
    '''
    Calculated the expected p-value to observe in potential SASEs
    
    requires
        long_fn: long output file name (contains potential SASEs)
        EVR_fn: evr file name
    yields
        local_expected: array of expected values
    '''
    
    long_generator = generate_long(long_fn)
    long_chrom, long_start, long_end, long_pvalue = long_generator.next()
    long_overlap = [defaultdict(int)]
    
    evr_lengths = defaultdict(int)
    evr_sums = defaultdict(int)
    previous_chrom = None
    seen = set()
    
    for line in open(EVR_fn, "r"):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        label = fields[3]
        seen.add(chrom)
        
        if previous_chrom == chrom or previous_chrom == None:
        
            evr_lengths[label] += end - start
        
            while chrom == long_chrom and long_start < end:
                if long_end >= start:
                    overlap = min(end, long_end) - max(start, long_start)
                    long_length = float(long_end - long_start)
                    evr_sums[label] += (long_pvalue / long_length) * overlap
                    long_overlap[-1][label] += overlap 
                try:
                    long_chrom, long_start, long_end, long_pvalue = long_generator.next()
                    if chrom == long_chrom:
                        long_overlap.append(defaultdict(int))
                except StopIteration:
                    break
                                        
        else:            
            local_expected = calc_expected(evr_lengths, evr_sums, long_overlap)
                        
            #print previous_chrom, local_expected
            yield local_expected
            
            long_overlap = [defaultdict(int)]
            evr_lengths = defaultdict(int)
            evr_sums = defaultdict(int)
            
            evr_lengths[label] += end - start
        
            while chrom == long_chrom and long_start < end:
                if long_end >= start:
                    overlap = min(end, long_end) - max(start, long_start)
                    long_length = float(long_end - long_start)
                    evr_sums[label] += (long_pvalue / long_length) * overlap
                    long_overlap[-1][label] += overlap 
                try:
                    long_chrom, long_start, long_end, long_pvalue = long_generator.next()
                    if chrom == long_chrom:
                        long_overlap.append(defaultdict(int))
                except StopIteration:
                    break
                    
        previous_chrom = chrom
    
    #print previous_chrom, local_expected 
    local_expected = calc_expected(evr_lengths, evr_sums, long_overlap)
    yield local_expected
        
