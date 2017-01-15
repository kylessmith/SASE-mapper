


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


def mut_intersect(muts_fn, evrs_fn):
    
    muts_generator = record_muts(muts_fn)
    mut_chrom, mut_pos, mut_sample = muts_generator.next()
    
    for line in open(EVR_fn, 'r'):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        label = fields[3] #make input
        
        while chrom == mut_chrom and mut_pos < end:
            if mut_pos >= start:
            
            try:    
                mut_chrom, mut_pos, mut_sample = muts_generator.next()
            except StopIteration:
                break                    