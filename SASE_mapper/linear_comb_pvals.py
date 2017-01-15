import numpy as np
from tpm import tpm
import bisect, time


def combine_pvals(sample_pos, sample_pvals, region_start, region_end, n_samples, tau=0.05):
    '''
    requires
        sample_starts: dictionary of numpy arrays of mutation positions
        sample_pvals: dictionary of numpy arrays of pvalues corresponding to sample_starts
        chrom_length: length of current chromosome
        chrom_start: position to start at (default = 0)
        tau: tau value for truncated product method
    returns
        chrom_c_pvals: numpy array of combined pvalues
        starts: numpy array of start positions corresponding to chrom_c_pvals
        ends: numpy array of end positions corresponding to chrom_c_pvals
    '''
    
    #create array the length of the region
    i = np.arange(region_start, region_end+1)
    #initialize array to hold pvalues
    total_pvalues = np.ones( (len(sample_pos.keys()), region_end-region_start) )
    #array to keep track of positions containing the same pvalues
    indices = np.zeros( region_end+1-region_start )
    
    chrom_c_pvals = []
    starts = []
    ends = []
    
    for k, sample in enumerate(sample_pos):
        sample_pvals[sample] = np.append(sample_pvals[sample], 1.0) #append 1.0 so -1 in pval_indices index to 1.0
        #find which base positions are assigned which pvalues
        pval_indices = np.searchsorted(sample_pos[sample], i, side='right')
        pvals = sample_pvals[sample][pval_indices-1]
        #assign pvalues to base positions
        total_pvalues[k] = pvals[:-1]
        indices += pval_indices
    
    #total_pvalues = total_pvalues.T
    previous_indices = np.nan
    
    for k in xrange(total_pvalues.shape[1]):
        #if indices are the same as the previous ones then do not calculate pvalue again
        if indices[k] != previous_indices:
            pval_list = np.zeros(n_samples)
            pval_list[:len(total_pvalues[:,k])] = total_pvalues[:,k]
            comb_pval = tpm(tau, 1000, time.time(), list(pval_list))
            chrom_c_pvals.append(comb_pval)
            starts.append(i[k])
            
        previous_indices = indices[k]
    
    starts = np.array(starts)
    ends = np.append(starts[1:], region_end)
        
    return np.array(chrom_c_pvals), starts, ends
        