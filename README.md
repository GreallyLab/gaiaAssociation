# Gaia Association

## Overview:

The gaiaAssociation package performs Regulatory LandscapeEnrichment Analysis (RLEA). RLEA tests for enrichment of sets
of loci in cell type-specific open chromatin regions (OCRs) in the genome. Given chromatin accessibility data (e.g. ATAC-seq, DNase-seq) and loci sets (e.g. de-novo mutations, rare variants, GWAS) gaiaAssociation will detect any cell-specific enrichment of loci.

A detailed guide and example notebook for using Gaia within a jupyter notebook environment can be found [here](https://github.com/GreallyLab/gaiaAssociation-Example-Guide).

![alt text](https://github.com/samrosean/images/blob/main/logo_with_border_transparent.png)

## Author/Support

Samuel Rosean, samuel.rosean@einsteinmed.edu --- https://github.com/samrosean

****

## Step 0: Installing Python and Pip:

To use gaiaAssociation as a command line tool or as a python library, the user will need python and pip installed onto their computer. To install python you can follow the instructions for your particular operating system given by Python:

 https://www.python.org/downloads/

Installing python from python.org should ensure that pip is already installed and ready to use. If this is not true, or you have installed python through another method you can follow this guide from pip’s documentation:

 https://pip.pypa.io/en/stable/installation/

The above link will also help if you need to upgrade or update your pip installation.

### Installing Required Libraries
		
gaiaAssociation requires seven python dependancies: scipy, pandas, numpy, pyranges, seaborn, matplotlib, and setuptools. When installing gaiaAssociation through pip it will attempt to resolve these dependencies and install the necceassary versions. If this fails, each package can be downloaded individually using pip, for example:

    pip install scipy
		
## Step 1: Install gaiaAssociation
	
### PyPi (Recommended)

#### For Command Line Usage:

To install gaiaAssociation from pypi, run a pip command to install it from the pypi distribution archive

	pip install gaiaAssociation

To get info on how to run from the command line run:

	gaia --help

 #### For Python Usage:
 
 A detailed guide and example notebook for using Gaia within a jupyter notebook environment can be found [here](https://github.com/GreallyLab/gaiaAssociation-Example-Guide).

 If you wish to use gaiaAssociation within a jupyter notebook or within a python script or notebook of your choice, then install using pip:

 	!pip install gaiaAssociation

Then, you will want to import the main function into your workspace

	import gaiaAssociation.gaiaAssociation

This will allow you to use gaia inline as a python function:

	gaiaAssociation.gaiaAssociation.gaiaAssociation("User/OCRfiles", "/User/lociFiles", "/User/chrsize.csv", "/User/Output", lociSelection = "/Users/lociGroups.tsv", windowSize = 10000)

### Github

You can also download the source package from this github repository, or through the terminal using the command given below (though we only recommend this method for people with experience building python packages and managing python environments):

	wget https://github.com/GreallyLab/gaiaAssociation/archive/main.zip

And then unpacking this zip file using the command.

	unzip main.zip

From the location of this newly installed copy of the gaiaAssociation repository, run the setup command

	python setup.py install

Gaia will now be runnable from the command line as described above within this folder.

## Step 2: Using gaiaAssociation:

A detailed guide and example notebook for using Gaia within a jupyter notebook environment can be found [here](https://github.com/GreallyLab/gaiaAssociation-Example-Guide).


When using gaiaAssociation on the command line you can check your installation and to get basic information on using gaiaAssociation by typing this command into your terminal.

	gaia --help

If correctly installed it should print out useful information about each variable and the flag to associate with each. gaiaAssocation has four required arguments, and five optional arguments. 

#### A Basic Run

A default run on the command line using just the required arguments will therefore look like:

```
gaia \
-a user/documents/atac \
-g user/documents/loci \
-c user/chrom/chrsize.csv \
-o user/documents/output
```

And a default run in a python notebook will look like:

```
gaiaAssociation.gaiaAssociation.gaiaAssociation("user/documents/atac", "user/documents/loci", "user/chrom/chrsize.csv", user/documents/output")
```

The four required arguments, which define a basic run are:

#### Required Arguments

OCR Folder: the folder location of the OCR bed files stored in .txt format. The first three columns should be labeled “Chromosome”, “Start”, and “End” in that order.

	-a, --atac  (either flag will work)

Loci Folder: the folder location of the loci files stored in .tsv and/or .csv format.
These files can be formatted differently depending on the type of loci being compared. If these loci are a single base pair long then only two columns are required: “CHR_POS” and “CHR_ID” which refer to its genome location and its chromosome name respectively. If these are variably sized loci you can include “Start” and “End” instead of “CHR_POS”. If you would like to give specific labels to each loci set, since by default we will name the loci set by their filename, you can include a column titled "DISEASE/TRAIT", which allows for multiple loci sets to be analyzed within one file.

	-g , --loci (either flag will work)

Chromosome Size: The location of a chromosome size file stored in a .csv format. This should only have two columns “Chromosome” and “size in bp”. This file can be used to subset what chromosomes you are interested in looking at enrichment within. If you only include chr1 and chr2, then enrichment will only be done relative to these two chromosomes. It is important to give it a path including the file, not just the path. (e.g. user/chrom/chrsize.csv not just user/chrom)

	-c, --chrom (either flag will work)

Output Folder: The folder location you want the results to be output into. If this folder does not already exist gaia will attempt to make it. If it does not have the permissions to do so it will exit and the user will have to run it with folder creating permissions, or they will have to make the folder themselves.

	-o, --output (either flag will work)

#### Optional Arguments

Peak Uniqueness: a cutoff value for OCR uniqueness (e.g. if given 12, then any atac peak found in more than 12 atac sets will be removed from all of them) - by default uniqueness is not considered).

```
  -u, --peakUniqueness
```

Loci Cutoff: a loci cutoff value, when given gaia will only consider loci groups (phenotypes or cohorts) with more loci than this cutoff value - by default this cutoff value is 0.

```
  -l, --lociCutoff
```

Specific Loci: a tsv or csv file with the specific loci groups you would like to use. This can be very helpful if using a large loci set with with many phenotypes, and you want to sort by more than just loci count. The tsv should have column names which also exist in the loci files, so that gaia can subset the loci based on these column values (e.g. if you have a column in your loci set titled "runs" and you would like to only use run 1, then a tsv with one column "runs" and one row value "1" will accomplish this). If the tsv/csv file includes multiple columns, say "STUDY", and "DISEASE/TRAIT", then the loci set will be subset specifically by both columns at once. So only those that match every value within a row will be considered as a group (e.g. if you subset by "runs" and by "patient", then if only those loci that match both of these values will be considered). (Note: this argument will supersede the loci cutoff argument, as all selected loci sets will be analyzed regardless of individual count)

```
  -s, --specificLoci
```

Masking Region: a bed file in a .txt format containing a set of regions that you want to subset every OCR region by. For example, a set of regions around the TSSs of a list of particular genes. This will reduce the OCR regions to just those that overlap with this given set of regions. This can be used to detect cell-specific + site-specific enrichment differences.

```
  -m, --maskRegion
```

Window Size: an integer given to represent the size of windows in bp that the user would like to divide the chromosome into. This method is based on the sinib tool (https://github.com/boxiangliu/sinib), which requires the chromosome be divided into a series of equal length windows to be able to model them as a series of binomials. The default value is 100,000 bp, but this value can be changed to increase specificity or decrease sepcificity. The function of the window size is to only consider the local environment when determing loci enrichment, so consideration should be made to what the user considers local in their particualr context.

```
  -w, --windowSize
```

A run from the command line using these flags will therfore look like:

```
gaia \
-a user/documents/atac \
-g user/documents/loci \
-c user/chrom/chrsize.csv \
-o user/documents/output \
-l 2000 \
-u 10 \
-m user/documents/mask.txt \
-w 1000000
```

Within a python file it will look like:

```
gaiaAsscoiation("user/documents/atac", "user/documents/loci", "user/chrom/chrsize.csv", "user/documents/output", lociCutoff = 2000, peakUniqueness = 10, maskRegion = "user/documents/mask.txt", windowSize = 1000000)
```
## Additional Notes

1. Gaia is designed to be run on tsv and csv files in a utf-8 encoding format. Problems may arise from using files encoded in a different format, we recommend you ensure files are saved in this encoder if errors are occuring.

2. Gaia does not have a built in genome build. Since the user provides the chromosome sizes, the OCR bed files, and the loci, as long as all three of these files are based on the same genome build it will run regardless of species/genome build. However, this does require that the user be vigilant, as results will still be generated if you compare OCRs from the a different genome build to the loci, these results will just be biologically meaningless.

3. Files inputted into Gaia do not need to be sorted. The loci files or bed files do not need to be sorted for Gaia to properly run.


## How Gaia Works:

By dividing each chromsome into roughly equivalent gaia window sizes based on a user-given value, enrichment is modeled as a binomial variable for each window where loci are found (wherein the probability is determined by the proportion of the window covered by open chromatin regions and the count is number of loci found within that window). The sum of these binomial variables are compared against the number of global overlaps between a cell-type's OCRs and the given loci set. The non-identical binomial variables are summed utilizing the method developed by Boxiang Liu and Thomas Quertermous (https://journal.r-project.org/archive/2018/RJ-2018-011/RJ-2018-011.pdf).

### Version 1.2.1
