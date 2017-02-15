from pybedtools import BedTool, Interval
from collections import defaultdict, OrderedDict
import numpy as np
import sys

def merge(overlaps, b, min_samples):
    '''
    requires
        overlaps: BedTool of intervals
        b: distance between overlaps to merge
        min_samples: minimum overlapping samples to report
    yields: tuple(current_chrom: current chromosome,
                  current_start: current start position,
                  current_end: current end position, 
                  number of samples, 
                  , separated string of samples)
    '''
    
    current_chrom = None
    current_start = 0
    current_end = 0
    samples = set()
    
    for overlap in overlaps:
        chrom = overlap[0]
        start = int(overlap[1])
        end = int(overlap[2])
        
        overlap_samples = set(overlap[4].split(","))
        
        if current_chrom == chrom or current_chrom == None:
            if (current_end + b) >= start:
                current_end = end
                samples = samples.union(overlap_samples)
            elif current_end != 0:
                if len(samples) >= min_samples:
                    yield current_chrom, current_start, current_end, len(samples), ','.join(list(samples))
                current_chrom = chrom
                current_start = start
                current_end = end
                samples = overlap_samples
            else:
                current_chrom = chrom
                current_start = start
                current_end = end
                
        else:
            if len(samples) >= min_samples:
                yield current_chrom, current_start, current_end, len(samples), ','.join(list(samples))
            current_chrom = chrom
            current_start = start
            current_end = end
            samples = overlap_samples
    
    if len(samples) >= min_samples:
        yield current_chrom, current_start, current_end, len(samples), ','.join(list(samples))


def generate_overlaps(overlaps, min_samples, b):
    '''
    create BedTool for overlapping regions
    
    requires
        overlaps: BedTool of intervals
        b: distance between overlaps to merge
        min_samples: minimum overlapping samples to report
    yields
        overlap: pybedtools.Interval object
    '''
    
    #if number of samples overlapping region is above min_samples,
    #yield interval
    for overlap in overlaps:
        
        if int(overlap[3]) >= min_samples:
            yield overlap


def cat_multi_beds(bed1, bed2, bed2_shift):
    '''
    concatenate 2 BedTools
    
    requires
        bed1: pybedtools.BedTool object
        bed2: pybedtools.BedTool object
        bed2_shift: 
    yield: string of tab delimited results
    '''
    
    fmt = "{chrom}\t{start}\t{end}\t{n_samples}\t{samples}"
    
    for fields in bed1:
        chrom = fields[0]
        start = fields[1]
        end = fields[2]
        n_samples = fields[3]
        samples = fields[4]
        yield fmt.format(**locals())
    
    for fields in bed2:
        chrom = fields[0]
        start = fields[1]
        end = fields[2]
        n_samples = fields[3]
        temp_samples = fields[4]
        
        samples = ','.join([str(int(sample)+bed2_shift) for sample in temp_samples.split(",")])
        yield fmt.format(**locals())
        

def filter_multi_bed(new_multi, min_samples):
    '''
    requires
        new_multi: pybedtools.BedTool object
        min_samples: minimum overlapping samples to report
    yields: tuple(chrom name, start, end, number of samples)
    '''
    
    #if number of samples overlapping region is above min_samples,
    #yield interval
    for region in new_multi:
        if region[3] >= min_samples:
            yield region[0], region[1], region[2], region[3]


def merge_multi_beds(multi_beds, b, min_samples):
    '''
    requires
        multi_beds: list of BedTools
        b: distance between overlaps to merge
        min_samples: minimum overlapping samples to report
    yields
        s: merged BedTool
    '''
    
    new_multi = []
    
    for i in xrange(len(multi_beds)-1):
        if i == 0:
            bed1 = multi_beds[i]
        else:
            bed1 = BedTool(new_multi)

        bed2 = multi_beds[i+1]
        bed2_shift = (i+1)*200
        
        #combine bed1 and bed2
        bed1_and_bed2 = BedTool(cat_multi_beds(bed1, bed2, bed2_shift)).saveas()
        if len(bed1_and_bed2) > 0:
            new_multi = merge(bed1_and_bed2.sort(), b, 1)
    #if number of samples overlapping region is above min_samples,
    #yield interval
    for s in new_multi:
        if s[3] >= min_samples:
            yield s


def many_multiIntersect(bedtool_names, min_samples, b):
    '''
    requires
        bedtool_names: list of BedTool names
        min_samples: minimum overlapping samples to report
        b: distance between overlaps to merge
    returns:
        overlaps: pybedtools.BedTool object
    '''
    
    n_bedtools = len(bedtool_names)
    multi_beds = []
    
    #conduct merges in chunks of 200
    for i in xrange(0,n_bedtools,200):
        try:
            bed_names = bedtool_names[i:i+200]
        except IndexError:
            bed_names = bedtool_names[i:]
        
        empty_bed = BedTool()
        multi_bed = merge(empty_bed.multi_intersect(i=bed_names), b, 1)
        multi_beds.append(multi_bed)
        
    overlaps = merge_multi_beds(multi_beds, b, min_samples)
    
    return overlaps


def identify_overlaps(pvalue_bedtools, min_samples=3, b=30):
    '''
    conduct multiIntersectBed to identify regions present in multiple samples
    return overlap region and sample pvalues
    
    requires
        pvalue_bedtools: list of BedTools
        min_samples: minimum overlapping samples to report
        b: distance between overlaps to merge
    returns
        overlaps: pybedtools.BedTool object
    '''
    
    bedtool_names = [bed.fn for bed in pvalue_bedtools.values()]
        
    #check if more than 200 BedTools will be produced,
    #if so run in chunks
    if len(pvalue_bedtools) > 200:
        overlaps = many_multiIntersect(bedtool_names, min_samples, b)
    
    else:
        empty_bed = BedTool()
        overlaps = merge(empty_bed.multi_intersect(i=bedtool_names), b, min_samples)
    
    return overlaps