from evr_permutation import evr_permute
import multiprocessing
import argparse
import time
from collections import OrderedDict
from find_overlaps import identify_overlaps
from linear_comb_pvals import combine_pvals
from tempfile import NamedTemporaryFile
from SFC import find_peaks
from pval_plot import plot
import functools
from pybedtools import BedTool


def combine_overlap_pvalues(overlap, tau, sample_muts, sample_pvals, n_samples):
    
    chrom_start = overlap[1]
    region_length = overlap[2]
    chrom_c_pvals, starts, ends = combine_pvals(sample_muts, sample_pvals, chrom_start, region_length, n_samples, tau)
        
    return starts, ends, chrom_c_pvals


def write_long(tau, sample_muts, sample_pvals, n_samples, out_long, p, overlap_generator, chrom):
    
    g = functools.partial(combine_overlap_pvalues, tau=tau, sample_muts=sample_muts, sample_pvals=sample_pvals, n_samples=n_samples)
    overlap_pvals = p.map(g, overlap_generator)
    
    for starts, ends, chrom_c_pvals in overlap_pvals:
        if len(starts) != 0:
            for i in xrange(len(chrom_c_pvals)):
                if chrom_c_pvals[i] < 1:
                    out_long.write(chrom+'\t'+str(starts[i])+'\t'+str(ends[i])+'\t'+str(chrom_c_pvals[i])+'\n')
                

def read_genome(genome_fn):
    
    genome = OrderedDict()
    
    for line in open(genome_fn, 'r'):
        fields = line.strip().split("\t")
        genome[fields[0]] = int(fields[1])
        
    return genome


def run_SASE_mapper(muts_fn, evrs_fn, genome_fn, cpu, tau, prefix, gap_bp,
                    min_pval, min_samples, max_intermut_pval, pth, local, global_sfc):
    
    out_short_fn = prefix+'_short.bed'
    out_long_fn = prefix+'_long.bed'
    temp_file = NamedTemporaryFile()
    #temp_file = open("pan_long_temp.bed", 'w')
    
    p = multiprocessing.Pool(cpu)
        
    #get dictionary {chrom:chrom_length,...}
    genome = read_genome(genome_fn)
    
    pvalue_bedtools_generator = evr_permute(muts_fn, evrs_fn, p, min_pval, max_intermut_pval)
    
    # think about while loop so you can delete as you go
    
    for chrom, n_samples, pvalue_bedtools, sample_pos, sample_pvals in pvalue_bedtools_generator:
        print chrom
        print "found", n_samples, "samples"
        if n_samples == 0:
            continue
        print "finding overlaps"
        overlap_generator = identify_overlaps(pvalue_bedtools, min_samples, gap_bp)
        print "combining pvalues"
        write_long(tau, sample_pos, sample_pvals, n_samples, temp_file, p, overlap_generator, chrom)
        print "done combining pvalues"
    
    temp_file.flush()
    print "writing long file"
    BedTool(temp_file.name).sort().saveas(out_long_fn)
    temp_file.close()
    
    print "finding peaks"
    find_peaks(out_long_fn, pth, out_short_fn, genome, muts_fn, min_samples, gap_bp, local, evrs_fn, global_sfc)
    print "plotting"
    plot(out_long_fn, genome, prefix)


def main():
    
    '''
    runs main ananlysis
    '''
    
    parser=argparse.ArgumentParser()
    
    parser.add_argument('--m', help='input mutation file', required=True)
    parser.add_argument('--e', help='input segmentation file', default='')
    parser.add_argument('--g', help='input genome file', required=True)
    parser.add_argument('--p', help='prefix for output files', required=True)
    parser.add_argument('--cpu', help='number of processors to use (default = max-1)', default = multiprocessing.cpu_count()-1, type=int)
    parser.add_argument('--dist_pval', help='maximum p-value for intermutation distances (default = 0.05)', default = 0.05, type=float)
    parser.add_argument('--tau', help='tau to use for combined p-value (default = 0.05)', default = 0.05, type=float)
    parser.add_argument('--b', help='maximum number of base pairs to combine regions (default = 30)', default = 30, type=int)
    parser.add_argument('--min_p', help='minimum detectable p-value (default = 1e-5)', default=1e-5, type=float)
    parser.add_argument('--min_samples', help='number of samples with mutations (default = 3)', default=3, type=int)
    parser.add_argument('--pth', help='p-value threshold for SFC, lower is stricter (default = 1e-6)', default=1e-6, type=float)
    parser.add_argument('--local', help='flag to run local SFC', action="store_true", default=False)
    parser.add_argument('--global_sfc', help='flag to run global SFC', action="store_true", default=False)
    
    args=parser.parse_args()
    
    
    genome_fn = args.g
    muts_fn = args.m
    evrs_fn = args.e
    cpu = args.cpu
    tau = args.tau
    prefix = args.p
    gap_bp = args.b
    min_pval = args.min_p
    min_samples = args.min_samples
    max_intermut_pval = args.dist_pval
    pth = args.pth
    local = args.local
    global_sfc = args.global_sfc
    
    start = time.clock()
    
    run_SASE_mapper(muts_fn, evrs_fn, genome_fn, cpu, tau, prefix, gap_bp,
                    min_pval, min_samples, max_intermut_pval, pth, local, global_sfc)
                        
    print time.clock() - start
    
        
if __name__ == "__main__":
    main()