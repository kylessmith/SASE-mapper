import numpy as np
import scipy.stats
from pybedtools import BedTool
from local_expected import calculate_local_expected
import sys


def generate_summary(short_temp, muts_fn, min_samples):
    '''
    requires
        short_temp: BedTool object
        muts_fn: mutation file name
        min_samples: minimum number of samples overlapping a region
    yields: tab delimited string of results
    '''
    
    fmt = "{chrom}\t{start}\t{end}\t{n_samples}\t{n_muts}\t{SFC}\n"
    muts_per_region = BedTool(muts_fn).intersect(short_temp, stream=True, wo=True)
    
    previous_region = None
    samples = set()
    n_muts = 0
    
    for i in muts_per_region:
        region = (i[-5], i[-4], i[-3], i[-2])
        if previous_region != region and previous_region != None:
            chrom = previous_region[0]
            start = previous_region[1]
            end = previous_region[2]
            SFC = previous_region[3]
            n_samples = len(samples)
            
            if n_samples >= min_samples:
                yield fmt.format(**locals())
            
            samples = set()
            n_muts = 0
        
        previous_region = region
        n_muts += 1
        samples.add(i[3])
    
    chrom = region[0]
    start = region[1]
    end = region[2]
    SFC = region[3]
    n_samples = len(samples)
    if n_samples >= min_samples:
        yield fmt.format(**locals())

def generate_short(chroms, starts, ends, sig_regions):
    '''
    requires
        chroms: array of chromosome names
        starts: array of start positions
        ends: array of end positions
        sig_regions: array of SFC scores
    yields: tab delimited string of results
    '''
    
    fmt = "{chrom}\t{start}\t{end}\t{SFC}\n"
    
    for i in xrange(len(starts)):
        SFC = sig_regions[i]
        if SFC > 0:        
            chrom = chroms[i]
            start = starts[i]
            end = ends[i]
            
            yield fmt.format(**locals())


def efc(n, X, expected, pth=1e-6):
    '''
    requires
        n: numpy array of lengths
        X: numpy array of signals
        expected: expected value
        pth: p-value to set intervals
    returns
        sig_regions: array of SFC scores
    '''
    ## Determine the upper and lower bounds
    ## using equation for Z in paper but solve
    ## for pe
    
    ## t = n* in paper
    #t = np.maximum(n,10)
    t = np.maximum(n, 1)
    Z = -scipy.stats.norm.ppf(pth)
    a = (2*n*X) + (t*Z**2)
    b = np.sqrt( (4*n*X) + (t*Z**2) - (4*X**2) )
    cp = np.sqrt(t)*Z
    cn = -np.sqrt(t)*Z
    d = 2*(n**2 + (t*Z**2))
    
    up_limit = (cp*b+a)/d
    low_limit = (cn*b+a)/d
    
    up_mean = n*up_limit
    low_mean = n*low_limit

    depleted = n*expected>up_mean
    enriched = n*expected<low_mean

    sig_regions = np.zeros(len(X))

    sig_regions[enriched] = np.log2(low_mean[enriched]/(n[enriched]*expected))
    sig_regions[depleted] = np.log2(up_mean[depleted]/(n[depleted]*expected))
        
    return sig_regions
    

def local_efc(n, X, expected, pth=1e-6):
    '''
    requires
        n: numpy array of lengths
        X: numpy array of signals
        expected: array of expected values
        pth: p-value to set intervals
    returns
        sig_regions: array of SFC scores
    '''
    ## Determine the upper and lower bounds
    ## using equation for Z in paper but solve
    ## for pe
    
    ## t = n* in paper
    #t = np.maximum(n,10)
    t = np.maximum(n, 10)
    Z = -scipy.stats.norm.ppf(pth)
    a = (2*n*X) + (t*Z**2)
    b = np.sqrt( (4*n*X) + (t*Z**2) - (4*X**2) )
    cp = np.sqrt(t)*Z
    cn = -np.sqrt(t)*Z
    d = 2*(n**2 + (t*Z**2))
    
    up_limit = (cp*b+a)/d
    low_limit = (cn*b+a)/d
    
    up_mean = n*up_limit
    low_mean = n*low_limit

    depleted = n*expected>up_mean
    enriched = n*expected<low_mean    

    sig_regions = np.zeros(len(X))

    sig_regions[enriched] = np.log2(low_mean[enriched]/(n[enriched]*expected[enriched]))
    sig_regions[depleted] = np.log2(up_mean[depleted]/(n[depleted]*expected[depleted]))
        
    return sig_regions


def parse_long(long_fn, pth):
    '''
    requires
        long_fn: long output file name
        pth: p-value to set SFC intervals
    returns
        array of chromosome names
        array of starting positions
        array of end positions
        array of pvalues
    '''
    
    chroms = []
    starts = []
    ends = []
    pvalues = []
    
    for line in open(long_fn, "r"):
        fields = line.strip().split("\t")
        pvalue = float(fields[3])
        
        if pvalue < 1:
            if pvalue < (pth*10):
                pvalue = (pth*10)
            chroms.append(fields[0])
            starts.append(int(fields[1]))
            ends.append(int(fields[2]))
            pvalues.append(-np.log10(pvalue))
        
    return np.array(chroms), np.array(starts), np.array(ends), np.array(pvalues)


def find_peaks(long_fn, pth, short_fn, genome, muts_fn, min_samples, b, seg_local, segs_fn, global_sfc):
    '''
    requires
        long_fn: long output file name
        pth: p-value to set SFC intervals
        short_fn: short output file name
        genome: dictionary of chromosome names -> chromosome length
        muts_fn: mutations file name
        min_samples: minimum number of samples to overlap a region
        b: distance to merge intervals
        seg_local: boolean to conduct local analysis
        segs_fn: evr file name
        global_sfc: boolean to conduct global SFC analysis
    writes results to file
    '''
    chroms, starts, ends, pvalues = parse_long(long_fn, pth)
    
    lengths = ends - starts
    coverage = lengths.sum()
    
    sig_regions = np.zeros(len(starts))
    if seg_local:
        expected_gen = calculate_local_expected(long_fn, segs_fn)
        
    if global_sfc:
        expected = np.sum(pvalues) / float(np.sum(genome.values()))
        
    for chrom in np.unique(chroms):
        print "\t", chrom
        if seg_local:
            while True:
                expected = expected_gen.next()
                if len(expected) == len(pvalues[chroms==chrom]):
                    break
        elif not global_sfc:
            expected = np.sum(pvalues[chroms==chrom]) / genome[chrom]
         
        if seg_local:
            sig_regions[chroms==chrom] = local_efc(lengths[chroms==chrom], pvalues[chroms==chrom], expected, pth)
        else:
            sig_regions[chroms==chrom] = efc(lengths[chroms==chrom], pvalues[chroms==chrom]/lengths[chroms==chrom], expected, pth)
        
    short_temp = BedTool(generate_short(chroms, starts, ends, sig_regions)).merge(d=b, c=4, o="mean")
    del lengths
    del chroms
    del starts
    del ends
    del pvalues
    short_out = BedTool(generate_summary(short_temp, muts_fn, min_samples)).saveas(short_fn)
    
    
    