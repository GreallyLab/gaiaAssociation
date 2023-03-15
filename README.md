# Gaia Association
## Version 0.1.0
## Samuel Rosean

****


## Overview:
Compare chromatin accessibility data (e.g. ATAC-seq, DNase-seq) against loci sets (e.g. de-novo mutations, rare variants) to detect cell-specific enrichment of loci.

## Guide:

### Installing Python and Pip:

To use gaiaAssociation as a command line tool, the user will need python and pip installed onto their computer. To install python you can follow the instructions for your particular operating system given by Python:

 https://www.python.org/downloads/

Installing python from python.org should ensure that pip is already installed and ready to use. If this is not true, or you have installed python through another method you can follow this guide from pip’s documentation:

 https://pip.pypa.io/en/stable/installation/

The above link will also help if you need to upgrade or update your pip installation.

### Installing Required Libraries
		
gaiaAssociation requires seven python dependancies: scipy, pandas, numpy, pyranges, seaborn, matplotlib, and setuptools. Each package can be downloaded using pip, for example:

    pip install scipy
		
## Install gaiaAssociation
	
		[Temporary]
### PyPi

To install gaiaAssociation for pypi, run a pip command to install it from the pypi distribution archive

	pip install -i https://test.pypi.org/simple/ gaiaAssociation==0.1.0
  
### Github

You can also download the source package from this github repository


From the location of this newly installed copy of the gaiaAssociation repository, run the setup command

 python setup.py install


## Using gaiaAssociation:

To check your installation and to get basic information on using gaiaAssociation type this command into your terminal.

	gaia --help

It should print out useful information about each variable and the flag to associate with each argument. gaiaAssocation has four required arguments, and six optional arguments. The four required arguments, which define a basic run are:

#### Required Arguments

ATAC Folder: the folder location of the ATAC bed files stored in .txt format. The first three columns should be labeled “Chromosome”, “Start”, and “End” in that order.

	-a, --atac  (either flag will work)

Loci Folder: the folder location of the loci files stored in .tsv format.
These files can be formatted differently depending on the type of loci being compared. If these loci are a single base pair long then only two columns are required: “CHR_POS” and “CHR_ID” which refer to its genome location and its chromosome name respectively. If these are variably sized loci you can include “Start” and “End” instead of “CHR_POS”. If you would like to give specific labels to each loci set, since by default we will name the loci set by their filename, you can include a column titled "DISEASE/TRAIT", which allows for multiple loci sets to be analyzed within one file.

	-g , --loci (either flag will work)

Chromosome Size: The location of a chromosome size file stored in a .csv format. This should only have two columns “Chromosome” and “size in bp”. This file can be used to subset what chromosomes you are interested in looking at enrichment within. If you only include chr1 and chr2, then enrichment will only be done relative to these two chromosomes. It is important to give it a path including the file, not just the path. (e.g. user/chrom/chrsize.csv not just user/chrom)

	-c, --chrom (either flag will work)

Output Folder: The folder location you want the results to be output into.

	-o, --output (either flag will work)

A default run will therefore look like:

```
gaia -a user/documents/atac -g user/documents/loci -c user/chrom/chrsize.csv -o user/documents/output
```

#### Optional Arguments

Peak Uniqueness: a cutoff value for ATAC uniqueness (e.g. if given 12, then any atac peak found in more than 12 atac sets will be removed from all of them) - by default uniqueness is not considered).

```
  -u, --peakUniqueness
```

Loci Cutoff: a loci cutoff value, when given gaia will only consider loci groups (phenotypes or cohorts) with more loci than this cutoff value - by default this cutoff value is 0.

```
  -l, --lociCutoff
```

Specific Loci: a txt file with the specific loci groups you would like to use. This can be very helpful if using a large loci set with with many phenotypes, and you want to sort by more than just loci count. (e.g. a complete database may include hundreds of loci groups, but you would like to manually select a subgroup of interest)

```
  -s, --specificLoci
```

Masking Region: a bed file in a a .txt format containing a set of regions that you want to subset each ATAC region by, for example a set of regions around the TSSs of particular genes of interest. This will reduce the ATAC regions to just those that overlap with this set of regions. This can be used to detect cell-specific + gene-specific enrichment.

```
  -m, --maskRegion
```

Print Loci: a flag which will save the overlapping LOCI for every single cell type for every single loci group as an independant file in a subfolder in the output, this is by default false

```
  -p
```

Z-Score Cutoff: a cutoff value for selecting which loci groups to keep after calculating the zscore enrichment for each loci for each given atac-seq set. Only loci groups with a zscore range (max-min) equal or greater to this cutoff will be kept. This can be useful for subsetting large groups of loci data into only those which have interesting insights, or unique cell enrichments.

```
  -z, --zcoreRangeCutoff
```

A run using these flags will therfore look like

```
gaia -a user/documents/atac -g user/documents/loci -c user/chrom/chrsize.csv -o user/documents/output -p -l 2000 -z 20 -u 10 -m user/documents/mask.txt
```
