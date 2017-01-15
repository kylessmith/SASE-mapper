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
    record the EVRs
    yield per chromosome
    '''
    
    for line in open(muts_fn, 'r'):
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        sample = fields[3]
        
        yield chrom, pos, sample


def record_muts_in_evr(muts_fn, evrs_fn):
    
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
    
    for line in open(evrs_fn, 'r'):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        label = fields[3] #make input
        seen.add(chrom)
        
        if chrom == previous_chrom or previous_chrom == None:
        
            try:
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            except KeyError:
                evr_lengths_cumsum[label] = [0]
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            evr_starts[label].append(start)
        
            while chrom == mut_chrom and mut_pos < end:
                if mut_pos >= start:
                    sample_mut_pos[mut_sample].append(mut_pos)
                    samples.add(mut_sample)
                    try:
                        evr_sample_muts[mut_sample][label] += 1
                    except KeyError:
                        evr_sample_muts[mut_sample][label] = 1
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
            while chrom != mut_chrom and mut_chrom in seen:
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
        else:
            yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, previous_chrom, list(samples)
            
            evr_lengths_cumsum = {}
            evr_starts = defaultdict(list)
            sample_mut_pos = defaultdict(list)
            evr_sample_muts = defaultdict(dict)
            samples = set()
            
            try:
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            except KeyError:
                evr_lengths_cumsum[label] = [0]
                evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end - start))
            evr_starts[label].append(start)
        
            while chrom == mut_chrom and mut_pos < end:
                if mut_pos >= start:
                    sample_mut_pos[mut_sample].append(mut_pos)
                    samples.add(mut_sample)
                    try:
                        evr_sample_muts[mut_sample][label] += 1
                    except KeyError:
                        evr_sample_muts[mut_sample][label] = 1
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
            while chrom != mut_chrom and mut_chrom in seen:
                try:    
                    mut_chrom, mut_pos, mut_sample = muts_generator.next()
                except StopIteration:
                    break
                    
        previous_chrom = chrom
                    
    yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, previous_chrom, list(samples)




def record_muts_in_evr2(muts_fn, evrs_fn):
    
    '''
    record intermut distance
    '''
    
    muts = BedTool(muts_fn)
    evrs = BedTool(evrs_fn)
    mut_evrs = evrs.intersect(muts, wao=True, sorted=True)
    
    evr_lengths_cumsum = {}
    evr_starts = defaultdict(list)
    sample_mut_pos = defaultdict(list)
    evr_sample_muts = defaultdict(dict)
    #evr_lengths = defaultdict(int)
    #evr_labels = []
    previous_chrom = None
    previous_evr_start = None
    samples = set()
    
    for i, evr in enumerate(mut_evrs):
        chrom1 = evr[0]
        start1 = int(evr[1])
        end1 = int(evr[2])
        label = evr[3]
        chrom2 = evr[4]
        start2 = int(evr[5])
        end2 = int(evr[6])
        sample = evr[7]
        
        if chrom1 == previous_chrom or previous_chrom == None:
        
            if start2 != -1:
                sample_mut_pos[sample].append(start2)
                samples.add(sample)
                try:
                    evr_sample_muts[sample][label] += 1
                except KeyError:
                    evr_sample_muts[sample][label] = 1
            
            if previous_evr_start != start1 or previous_evr_start == None:
                try:
                    evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end1 - start1))
                except KeyError:
                    evr_lengths_cumsum[label] = [0]
                    evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end1 - start1))
                evr_starts[label].append(start1)
        
        else:
            yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, previous_chrom, list(samples)
            evr_lengths_cumsum ={}
            evr_starts = defaultdict(list)
            sample_mut_pos = defaultdict(list)
            evr_sample_muts = defaultdict(dict)
            samples = set()
            
            if start2 != -1:
                sample_mut_pos[sample].append(start2)
                samples.add(sample)
                try:
                    evr_sample_muts[sample][label] += 1
                except KeyError:
                    evr_sample_muts[sample][label] = 1
            
            if previous_evr_start != start1 or previous_evr_start == None:
                try:
                    evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end1 - start1))
                except KeyError:
                    evr_lengths_cumsum[label] = [0]
                    evr_lengths_cumsum[label].append(evr_lengths_cumsum[label][-1] + (end1 - start1))
                evr_starts[label].append(start1)
                
        previous_evr_start = start1
        previous_chrom = chrom1
        
    yield evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, chrom1, list(samples)
    

def permute_mutations(sample, evr_lengths_cumsum, sample_muts, evr_starts, min_pval):
    
    total_muts =  np.sum(sample_muts[sample].values())
    n_perms = np.ceil((1/min_pval) / total_muts)+1
    shift = 0
    permutations = np.zeros((n_perms, total_muts))

    for label in sample_muts[sample]:
        label_mutations = sample_muts[sample][label]
        total_length = evr_lengths_cumsum[label][-1]
    
        permuted_muts = permutation(total_length, label_mutations, n_perms, np.array(evr_starts[label]), np.array(evr_lengths_cumsum[label]))
    
        permutations[:, shift:shift+permuted_muts.shape[1]] = permuted_muts
        shift += permuted_muts.shape[1]

    return np.ravel(np.diff(np.sort(permutations)))
    
    
def permutation(total_length, label_mutations, n_perms, evr_starts, evr_lengths_cumsum):
    
    permuted_muts = rand.random_int(total_length, (n_perms, label_mutations))
    
    index = evr_lengths_cumsum.searchsorted(permuted_muts, side="right") - 1
    permuted_muts -= evr_lengths_cumsum[index]
    permuted_muts += evr_starts[index]
    
    return permuted_muts
    

def get_mut_pvals(mut_positions, p, distribution=None, min_pval=1e-5,):
    
    '''
    requires
        mut_positions: numpy array of mutation positions
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
    
     return np.sum(y<=x)
     
     
def generate_bedtool(chrom, mut_pos, inter_mut_pvalues, max_intermut_pval):
    
    fmt = "{chrom}\t{start}\t{end}\t{pvalue}\n"
    
    for i in xrange(len(inter_mut_pvalues)-1):
        start = mut_pos[i]
        end = mut_pos[i+1]
        pvalue = inter_mut_pvalues[i]
        if pvalue <= max_intermut_pval:
            yield fmt.format(**locals())


def evr_permute(muts_fn, evrs_fn, p, min_pval, max_intermut_pval):
    muts_evr_generator = record_muts_in_evr(muts_fn, evrs_fn)
    
    # think about while loop so you can delete as you go
    
    for evr_lengths_cumsum, evr_starts, sample_mut_pos, evr_sample_muts, chrom, samples in muts_evr_generator:
        print "starting permutations for chrom", chrom
        #print [len(sample_mut_pos[s]) for s in sample_mut_pos]
        g = functools.partial(permute_mutations, evr_lengths_cumsum=evr_lengths_cumsum, sample_muts=evr_sample_muts, evr_starts=evr_starts, min_pval=min_pval)
        perm_distributions = p.map(g, samples)
        #print evr_sample_muts
        #perm_distributions = {sample:permute_mutations(sample, evr_lengths_cumsum, evr_sample_muts, evr_starts, min_pval) for sample in samples}
        
        sample_bedtools = {}
        n_samples = len(samples)
        sample_pvals = {}
        
        for i, sample in enumerate(samples):
            inter_mut_pvalues = get_mut_pvals(sample_mut_pos[sample], p, perm_distributions[i], min_pval)
            pval_bedtool = BedTool(generate_bedtool(chrom, sample_mut_pos[sample], inter_mut_pvalues, max_intermut_pval)).saveas()
            sample_bedtools[sample] = pval_bedtool
            sample_pvals[sample] = inter_mut_pvalues
            
        del perm_distributions
            
        yield chrom, n_samples, sample_bedtools, sample_mut_pos, sample_pvals