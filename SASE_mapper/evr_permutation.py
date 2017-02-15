from pybedtools import BedTool
import numpy as np
import rand
from collections import defaultdict
import multiprocessing
import functools
import warnings
import sys

warnings.filterwarnings('ignore')



def record_muts(muts_fn):
    '''
    muts_fn: file name of mutations BED file
    
    yields: tuple(chromosome name, position, sample name)
    '''
    
    #iterate through lines of mutations file
    for line in open(muts_fn, 'r'):
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        sample = fields[3]
        
        yield chrom, pos, sample


def record_muts_in_evr(muts_fn, evrs_fn):
    '''
    muts_fn: file name of mutations BED file
    evrs_fn: file name of equi-variant regions file
    
    yields: tuple(evr_lengths_cumsum: dictionary of arrays of length cummulative sums
                  evr_starts: dictionary of arrays of evr starting positions,
                  sample_mut_pos: dictionary of arrays of mutation positions,
                  evr_sample_muts: dictionary of dictionaries of the number of mutations in each evr,
                  previous_chrom: the  previous chromosome name,
                  samples: list of observed samples
    '''
    
    #initiate variables
    evr_lengths_cumsum = {}
    evr_starts = defaultdict(list)
    sample_mut_pos = defaultdict(list)
    evr_sample_muts = defaultdict(dict)
    previous_chrom = None
    samples = set()
    
    muts_generator = record_muts(muts_fn)
    mut_chrom, mut_pos, mut_sample = muts_generator.next()
    samples.add(mut_sample)
    seen = set()
    
    #interate through line of evr file
    for line in open(evrs_fn, 'r'):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        label = fields[3] #make input
        seen.add(chrom)
        
        #check if evr chromosome equals previous chromosome,
        #then add information on evrs
        if chrom == previous_chrom or previous_chrom == None:
            
            #check if evr label has been seen before
            try:
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            except KeyError:
                evr_lengths_cumsum[label] = [0]
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            evr_starts[label].append(start)
            
            #check if chromosome in evr match the mutation file chromosome name,
            #check if mutation position is less than evr end position
            while chrom == mut_chrom and mut_pos < end:
                #if mutation position is greater than the evr start positon,
                #then the mutation must be in the evr
                if mut_pos >= start:
                    sample_mut_pos[mut_sample].append(mut_pos)
                    samples.add(mut_sample)
                    try:
                        evr_sample_muts[mut_sample][label] += 1
                    except KeyError:
                        evr_sample_muts[mut_sample][label] = 1
                #interate through mutation generator,
                #if it is the end, then break
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
            
            #keep iterating through mutations until
            #mutation chromosome matches evr chromosome
            while chrom != mut_chrom and mut_chrom in seen:
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
        else:
            #change in chromosome detected, yield information
            yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, previous_chrom, list(samples)
            
            #reset variables
            evr_lengths_cumsum = {}
            evr_starts = defaultdict(list)
            sample_mut_pos = defaultdict(list)
            evr_sample_muts = defaultdict(dict)
            samples = set()
            
            #check if evr label has been seen before
            try:
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            except KeyError:
                evr_lengths_cumsum[label] = [0]
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            evr_starts[label].append(start)
            
            #check if chromosome in evr match the mutation file chromosome name,
            #check if mutation position is less than evr end position
            while chrom == mut_chrom and mut_pos < end:
                #if mutation position is greater than the evr start positon,
                #then the mutation must be in the evr
                if mut_pos >= start:
                    sample_mut_pos[mut_sample].append(mut_pos)
                    samples.add(mut_sample)
                    try:
                        evr_sample_muts[mut_sample][label] += 1
                    except KeyError:
                        evr_sample_muts[mut_sample][label] = 1
                #interate through mutation generator,
                #if it is the end, then break
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
            
            #keep iterating through mutations until
            #mutation chromosome matches evr chromosome
            while chrom != mut_chrom and mut_chrom in seen:
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
        previous_chrom = chrom
                    
    yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, previous_chrom, list(samples)


def permute_mutations(sample, evr_lengths_cumsum, sample_muts, evr_starts, min_pval):
    '''
    requires
        sample: sample name
        evr_lengths_cumsum: dictionary of arrays of length cummulative sums
        sample_muts: dictionary of dictionaries of the number of mutations in each evr
        evr_starts: dictionary of arrays of evr starting positions
        min_pval: minimum observable p-value from permutation
    returns: array of permuted mutation positions
    '''
    
    total_muts =  np.sum(sample_muts[sample].values())
    n_perms = np.ceil((1/min_pval) / total_muts)+1
    shift = 0
    permutations = np.zeros((n_perms, total_muts))
    
    #iterate through evr labels
    for label in sample_muts[sample]:
        label_mutations = sample_muts[sample][label]
        total_length = evr_lengths_cumsum[label][-1]
    
        permuted_muts = permutation(total_length, label_mutations, n_perms, np.array(evr_starts[label]), np.array(evr_lengths_cumsum[label]))
    
        permutations[:, shift:shift+permuted_muts.shape[1]] = permuted_muts
        shift += permuted_muts.shape[1]

    return np.ravel(np.diff(np.sort(permutations)))
    
    
def permutation(total_length, label_mutations, n_perms, evr_starts, evr_lengths_cumsum):
    '''
    requires
        total_length: sum of evr lengths
        label_mutations: number of mutations occurring in given evr
        n_perms: number of permutations to conduct
        evr_starts: array of evr start positions
        evr_lengths_cumsum: array of cummulative sums of evr lengths
    returns
        permuted_muts: 2D array of permutated mutation positions
    '''
    
    permuted_muts = rand.random_int(total_length, (n_perms, label_mutations))
    
    index = evr_lengths_cumsum.searchsorted(permuted_muts, side="right") - 1
    permuted_muts -= evr_lengths_cumsum[index]
    permuted_muts += evr_starts[index]
    
    return permuted_muts
    

def get_mut_pvals(mut_positions, p, distribution=None, min_pval=1e-5,):
    
    '''
    requires
        mut_positions: numpy array of mutation positions
        p: pool of processors from multiprocessing.Pool
        distribution: array of permutation intermutation distances
        min_pval: minimum observable p-value from permutation
    returns
        inter_mut_pvalues: numpy array of pvalues
    '''
    
    #record distance between mutations
    inter_mut_distances = np.diff(mut_positions)
    
    if distribution == None:
        #TODO
        pass
    
    else:
        g = functools.partial(compare, y=distribution)
        n_less = p.map(g, inter_mut_distances)
        inter_mut_pvalues = np.array(n_less)/float(len(distribution))
        inter_mut_pvalues[inter_mut_pvalues==0] = min_pval*0.1
    
    return inter_mut_pvalues
    

def compare(x,y):
    '''
    requires
        x: array of floats
        y: array of floats
    
    returns: total number of instances y <= x
    '''
    return np.sum(y<=x)
     
     
def generate_bedtool(chrom, mut_pos, inter_mut_pvalues, max_intermut_pval):
    '''
    requires
        chrom: chromosome name
        mut_pos: array of mutation positions
        inter_mut_pvalues: array of intermutation p-values
        max_intermut_pval: maximum p-value to include in analysis
    
    yields: string of tab delimited results
    '''
    
    fmt = "{chrom}\t{start}\t{end}\t{pvalue}\n"
    
    for i in xrange(len(inter_mut_pvalues)-1):
        start = mut_pos[i]
        end = mut_pos[i+1]
        pvalue = inter_mut_pvalues[i]
        if pvalue <= max_intermut_pval:
            yield fmt.format(**locals())


def evr_permute(muts_fn, evrs_fn, p, min_pval, max_intermut_pval):
    '''
    muts_fn: file name of mutations BED file
    evrs_fn: file name of equi-variant regions file
    p: pool of processors from multiprocessing.Pool
    min_pval: minimum observable p-value from permutation
    max_intermut_pval: maximum p-value to include in analysis
    
    yields: tuple(chrom: chromosome name,
                  n_samples: number of samples,
                  sample_bedtools: dictionary of samples -> pybedtools.BedTool objects
                  sample_mut_pos: dictionary of samples -> arrays of mutation positions
                  sample_pvals: dictionary of samples -> arrays of intermutation p-values)
    '''
    
    muts_evr_generator = record_muts_in_evr(muts_fn, evrs_fn)
    
    #iterate through chromosomes
    for evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, chrom, samples in muts_evr_generator:
        print "starting permutations for chrom", chrom
        #generate permutation distributions for each sample
        g = functools.partial(permute_mutations, evr_lengths_cumsum=evr_lengths_cumsum,
                              sample_muts=evr_sample_muts, evr_starts=evr_starts, min_pval=min_pval)
        perm_distributions = p.map(g, samples)
        
        sample_bedtools = {}
        n_samples = len(samples)
        sample_pvals = {}
        
        #create BedTool intermutation p-values for each sample
        for i, sample in enumerate(samples):
            inter_mut_pvalues = get_mut_pvals(sample_mut_pos[sample], p, perm_distributions[i], min_pval)
            pval_bedtool = BedTool(generate_bedtool(chrom, sample_mut_pos[sample], inter_mut_pvalues, max_intermut_pval)).saveas()
            sample_bedtools[sample] = pval_bedtool
            sample_pvals[sample] = inter_mut_pvalues
            
        del perm_distributions
            
        yield chrom, n_samples, sample_bedtools, sample_mut_pos, sample_pvals