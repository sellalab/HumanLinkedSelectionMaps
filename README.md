## Overview
This README file explains the basics of how to set up and run the `HumanLinkedSelectionMaps` code and replicate the results presented in our paper (https://www.biorxiv.org/content/10.1101/2021.07.02.450762v1). The code can be adapted to work with data from different species or using different human datasets not explored in our paper. For completeness, we present the full pipeline including the steps for processing data (annotations, genetic maps, etc.). However, if you are planning to work with a different set of annotations or species data, many of these steps can be skipped. Note that you will still need to format your input data for compatability with data processing steps. More detailed documentaion can also be found in the relevant python and bash files used to run the pipeline.

Everything was run using Python 2.7, so some compatibility changes may be needed for Python 3.x+. 

The basic steps are: 
  1. Using B-maps
  2. Downloading and installing all software
  3. Updating file paths in python and bash scripts
  4. Downloading public datasets used in analysis
  5. Formatting and pre-processing data
  6. Running the inference
  7. Analyzing the inference results
  <br/><br/>

## 1. Using B-maps
If you are just interested in using the B-maps that are referenced in our paper, these can be found in the `/Bmaps/` directory, along with relevant README files. These are formatted using the hg19 genome build.
<br/><br/>

## 2. Downloading and installing all software
First, download the entire contents of this repository, including all python code, shell scripts and the modified version of the source code for  McVicker's calc_bkgd program (McVicker et al., 2009) in `/debug_bkgd/`. You will probably need to install some dependencies to run calc_bkgd, depending on your system. This step can bit a bit of a pain, depending on your programming skill level. I have included a list of the required C libraries and some notes on installation that were helpful in my case in the `/debug_bkgd/notes-on-dependencies.txt` file.

Additional software that are used for data processing are `phastCons` and `phyloFit` from the PHAST suite of tools (Siepel and Haussler 2004; Siepel, Bejerano et al. 2005). 
<br/><br/>

## 3. Updating file paths in python and bash scripts
Both python and bash scripts have a number of embedded paths that are used to call either other scripts or for loading data. The vast majority of these in python are linked to the `root_dir` variable in the `/classes/runstruct.py` file, so updating this path to fit your own directory structure will carry through 99% of downstream paths. There is a small subset of paths to `/ifs/scratch/c2b2/gs_lab/dam2214/scratch` that you will need to search+replace to match your directory structure. 

Similarly, the bash scripts stored in the `/sh/` folder generally include paths to the python binary file from my system `py=/path/` and calls to the python program that will be called by the shell script `prog=/path/` that will need to be updated to match the current file structure.
<br/><br/>

## 4. Downloading public datasets used in analysis
All of the annotations, genetic maps, polymorphism data, etc., are downloaded from publicly available datasets (see Supplementary Material Section 2 for details https://www.biorxiv.org/content/10.1101/2021.07.02.450762v1). Some of these, like the 100-vertebrate alignment MAF files used to generate phastCons scores and the 1000 Genomes Project VCF files, are quite large, so be sure you have ~ 1TB available if you plan to download everything at once. Many of the larger files (like the 100-vertebrate alignment MAF files) will only be used once to create smaller processed datasets, so the large originals can be deleted. Alternatively, some of the processing can be done by streaming remote files into a pipe using `curl`, although you will need to a stable connection for this option (see, for example, `/1000_genomes/sh/recodesubpop.sh` script, which processes a stream of 1000 Genomes Project VCF file data directly).

Most of the file paths that are used for processing data (particularly all of the processed data used in the inference and results analysis) are found in the class `classes/filestruct.py`. You should duplicate the directory structure in this class on your own system (or modify it, depending on your preference) so that all of the files will be read/written in the correct location. Some file paths are also embedded in the relevant bash scripts or python files directly, so these must be updated accordingly. 

In general, searching `/Users/davidmurphy/`, `/ifs/data/` and `/ifs/scratch/` throughout the code can be performed to systematically update the directory structure.
<br/><br/>

## 5. Formatting and pre-processing data

