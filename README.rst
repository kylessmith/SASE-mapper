Software to identify regions of interest with a higher than expected number of mutations in
multiple samples. 

Data Format
===========

Mutations input files must be in BED format(first columns are chrome, start, stop, sample)::

    chr1    11873   14409	sample1

Simple Example:

    https://raw.githubusercontent.com/kylessmith/SASE-mapper/master/example/mutations.bed


NOTE: Most input files are assumed to be *sorted* BED files

Invocation
==========

Running the following command will result in a more detailed help message::

    $ python -m SASE_mapper -h

Gives::

	  --m M                 input mutation file
	  --e E                 input segmentation file
	  --g G                 input genome file
	  --p P                 prefix for output files
	  --cpu CPU             number of processors to use (default = max-1)
	  --dist_pval DIST_PVAL
	                        maximum p-value for intermutation distances (default =
	                        0.05)
	  --tau TAU             tau to use for combined p-value (default = 0.05)
	  --b B                 maximum number of base pairs to combine regions
	                        (default = 30)
	  --min_p MIN_P         minimum detectable p-value (default = 1e-5)
	  --min_samples MIN_SAMPLES
	                        number of samples with mutations (default = 3)
	  --pth PTH             p-value threshold for SFC, lower is stricter (default
	                        = 1e-6)
	  --local               flag to run local SFC
	  --global_sfc          flag to run global SFC

QuickStart
==========

If your files are in sorted BED format, want to run on default parameters,
and want SFC calculated locally.


mutations in chromosome 21 for melanoma
---------------------------------------
::

    $ python -m SASE_mapper \
        --m example/mutations.bed \
        --e example/EVRs.bed \
        --g example/hg19.genome \
		--p example/test \
        --local \

The output will be shown in the following files::

	examples/test_short.bed
	examples/test_long.bed
	examples/test.png

Installation
============

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  

    https://continuum.io/downloads

Dependencies can be installed by::

    pip install -r requirements.txt

SASE_mapper also depends on BEDTools which is available from https://github.com/arq5x/bedtools2/
