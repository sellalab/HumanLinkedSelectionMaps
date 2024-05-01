## Overview
This README file explains the basics of how to set up and run the `HumanLinkedSelectionMaps` code and replicate the results presented in [our paper](https://elifesciences.org/articles/76065). The code can be adapted to work with data from different species or different human datasets not explored in our paper (see "Advanced usage" section). If you are just interested in downloading and using the B-maps, these can be found in the `Bmaps/` directory, along with relevant README files. B-maps are formatted using the hg19 genome build.

The minimal data required to estimate the effects of linked selection in a given population are:
- polymorphism data
- a list of neutral sites
- a genetic map (ideally not based on polymorphism data, which can be confounded with effects of linked selection)
- estimates of the neutral substitution rate with fairly high resolution (e.g., 6kb windows)
- annotations of candidate targets of purifying selection (optional)
- annotations of candidate targets of selective sweeps (optional)

Note that while either selection target annotation is optional, at least one of them is required: either targets of purifying selection resulting in background selection AND/OR targets of selective sweeps. See Appendix 1 Sections 2 and 3 for details on how we constructed the dataset for the results presented in our paper. In the sections below, we detail the directories where each type of file can be found in our downloadable dataset.

The following "Quick start" guide will allow you to repdroduce the inference results presented in the paper (i.e. the CADD 6% and phastCons 6% results shown in the figures).

**NOTE: The code in this repository was run using Python 2.7, so some compatibility changes may be needed for Python 3.x+.**


# QUICK START GUIDE
1. Download software and update `root_dir`
2. Download data
3. Create new `RunStruct` file
4. Compress data
5. Run inference
6. Create figures

## 1. Download software and update `root_dir`
First, download the entire contents of this repository. Name the folder `lsm_python` (or whatever you like). Paths to input and output files are set in the `RunStruct` class, which indlues a hard-coded `root_dir` variable at the top of the file `lsm_python/classes/runstruct.py`. This was set to have two options, a local option and a remote option, depending on the User path on your system. Both paths must be updated here:

```
# ROOT DIRECTORIES FOR MY COMPUTER AND CLUSTER
if os.getcwd().startswith('/Users/davidmurphy'):
    root_dir = '{}/linked_selection'.format(goog_dir)
else:
    root_dir = '/ifs/data/c2b2/gs_lab/dam2214/linked_selection'
```

**NOTE: Don't forget to add the full path to `lsm_python` to your `PYTHONPATH` variable in the terminal. This will enable modules to be imported within python code.**

## 2. Download data
Download the following [google drive folder](https://drive.google.com/drive/folders/1YJn54EXnPRwd1rp8-j-CpknAglUu5gE5?usp=drive_link), which includes all of the data to recreate the main results. To preserve the expected file structure, the contents of the folder must be put into the same directory where you've put `lsm_python`. This folder contains files processed to the second to last step required before running the inference. The final step, compressing the data into "bins", is run by the user (see Appendix 1 Section 1.3 for description and rationale for binning). 

## 3. Create new `RunStruct` file
The `RunStruct` class includes a representation of the complete file structure used for input and output as well as inference settings, inference results, and other important variables. The code for this class is found in `lsm_python/classes/runstruct.py` (see above section on the root directory). For a given configuration of the inference, a new `RunStruct` is initialized and saved to a text file that is used as the input to a number of other programs. The file `lsm_python/likelihood/initialize_runstruct.py` is set up to initialize a `RunStruct` and save it in text format for the inference results presented as "CADD 6%" in the main paper. It will also create some folders for initialization and results files. Run this code and check that the file `/<ROOT_DIR>/result/init_files/YRI.cadd94_gmask_v1.6_without_bstat.BS1.6.CS0.0.NOT_STARTED.initial.txt` was created in your root directory. 

To initialize a `RunStruct` for the other main result shown in the paper (phastCons 6%), change the line `b_an_1 = 'cadd94_gmask_v1.6_without_bstat'` to `b_an_1 = 'fish_cons94_new'` in `lsm_python/likelihood/initialize_runstruct.py` and run the code again. This should create a comparable text file with the new annotation in the file name.

The `RunStruct` is a multilayered data structure which includes several sub-classes as well, including the `FileStruct` which contains the full file structure for input and output files and the `OptimizerParams` which contains the initial and final parameters of the optimization.  

## 4. Compress data
As detailed in Appendix 1 Section 1.3, we compress neutral sites into "bins" where the values of exogenous parameters or pre-calculated "raw" B values are constant, which vastly downscales the number of calculations needed in the inference without any loss of precision in our predictions. The file `lsm_python/precalc/lh_inputs.py` will perform this step by taking as input the `RunStruct` file created in the previous step. The processing pipeline takes one chromosome at a time and can be used with external command line arguments to run each autosome in parallel on a cluster (see code segments in the `main` function that are commented out). In the current configuration of the `main` function, a `for` loop will just go through each chromosome in serial. Note that this step may take 10-20 minutes (depending on your system).

This step builds the entire data structure needed to run the inference and includes the key files listed here:
- Zipped numpy arrays of neutral sites (neutral "mask" `.npz` files) for each chromosome. These arrays are the length of the chromosome in their filename, with a 1 for neutral sites and a 0 for non-neutral sites (see more details of the process of curating neutral sites in Appendix 1 Sections 2 and 3 in our paper). These files are found in the folder `data/nmsk` in the Google Drive folder linked above.
- Pre-calculated maps of the effects of background selection for each chromosome and annotation (optional). These text files are found in the `precalc/cadd94_gmask_v1.6_without_bstat` directory for the current example. To calculate these files directly from genetic maps and annotated targets of background selection, see below "ADVANCED USAGE: 3. Pre-calculating background selection maps using `calc_bkgd`".
- Pre-calculated maps of the effects of selective sweeps for each chromosome and annotation (optional). These text files will be found in the `precalc/` directory but are not used in the current example. To calculate these files directly from genetic maps and annotated targets of sweeps, see below "ADVANCED USAGE: 4. Pre-calculating selective sweep maps".
- Estimates of the neutral substitution rate in bins of 6000 neutral sites. These are used to scale predictions of the effects of linked selection on diveristy levels by accounting for variation in the neutral substitution rate along the genome (see Appendix 1 Section 3 of our paper for details). These text files are found in the `data/phast/aligned_8_win_6000` directory for the current example.
- Polymorphism data from the population of interest (YRI in the current example). These are tables of counts of reference and alternate alleles at all biallelic SNPs in the 1000 Genomes Project Phase 3 dataset and are found in the `data/snps/frq/YRI` directory for the current example. This polymorphism data will be filtered in combination with the neutral mask files described above to only include neutral polymorphic sites. 

## 5. Run inference
As described in Appendix 1 Section 1.4, we run the inference in two optimization steps. In the first step, we run the optimization from 15 different starting conditions and keep the 3 results with the greatest likelihood. In the second step, we take the mean parameters of these 3 "best" runs from stage 1 and use these as the input parameters to the optimization. 

To run the first step, run the script `lsm_python/likelihood/inf_step_1.py`. This will read the `RunStruct` file (created above) and loop through the 15 different starting conditions (indexed 0-14), leading to 15 different results files in the folder `<ROOT_DIR>/result/final_files/cadd94_gmask_v1.6_without_bstat/`, with index 0-14 in the file labels. Note that this step is much faster if you can run each of the 15 optimizations in parallel and take the index as an external command line argument (see code segments in the `main` function that are commented out).

After all 15 optimizations finish, step 1 is complete and you can move on to step 2. To do so, run the script `lsm_python/likelihood/inf_step_2.py`. As in the above steps, this script is hard coded to continue running through the example for CADD 6% which we've started in the previous steps. The `main` function will just read the name of the results folder (i.e., `cadd94_gmask_v1.6_without_bstat`) as input. After the step 2 optimization is complete, a new file with the suffix `.composite.txt` will be saved to the folder `<ROOT_DIR>/result/final_files/cadd94_gmask_v1.6_without_bstat/`. This file contains the best-fitting parameters the current model configuration and data.

## 6. Create figures
After the inference is run and the final `.composite.txt` file has been created, you can use the inferred parameters to perform the analyses presented in the paper. Here we will cover the figures from the main paper (Figs. 2-5); to generate figures from the Appendix, see **ADVANCED USAGE** below.

Each of the figures requires an analysis step that results in a small table of values that are used to render the plots in Figs. 2-5. To simplify, each of these steps has been put into a single python file: `lsm_python/figures/main_figs_analysis.py`. Running this script will run each analysis for Figs. 2-5 and save the results in text or numpy array format to the folder `cadd94_gmask_v1.6_without_bstat`. You can explore each of the called analysis functions in more detail if you'd like to better understand the analyses, or change certain paramters (e.g., you could change the chromosome in figure 2a from chr1 to chr5).

To plot the figures, use the script `lsm_python/figures/main_figs_quickstart.py`. This script is set up like a python notebook with individual cells for each figure. The individual cells will make each of the plots (Figs. 2-5). You can also just run the whole script and generate all of the figures. Saved files will be written to `<ROOT_DIR>/result/final_files/mainfigs`.

**NOTE: The `results/final_files/` folder includes some additional results folders that are used for mutation rate estimates plotted in figure 4.**


## ADVANCED USAGE
1. Download and install third-party software used for data processing
2. Download publicly available datasets used for the inference 
3. Pre-calculating background selection maps using `calc_bkgd`
4. Pre-calculating selective sweep maps

## 1. Download and install third-party software used for data processing
To run the inference starting from scratch, you will need to generate new "raw" maps of the effects of background selection for a grid of selection effects using the modified version Graham McVicker's `calc_bkgd` program (McVicker et al., 2009) found in `debug_bkgd/`. To compile this C program, you will most likely need to install some dependencies, which may vary depending on your system, and this step can bit a bit of a pain, depending on your programming skill level. I have included a list of the required C libraries that were needed to get this compiled on MacOS and some notes on installation that were helpful in my case in the `/debug_bkgd/notes-on-dependencies.txt` file.

To identify conserved regions and also estimate substitution rates, data were processed with both the `phastCons` and `phyloFit` tools from the PHAST suite (Siepel and Haussler 2004; Siepel, Bejerano et al. 2005) found [here](http://compgen.cshl.edu/phast/index.php). These programs are well-documented very easy to install and use.

To process raw VCF files from the 1000 Genomes Project, we used [VCFTools](https://vcftools.sourceforge.net/man_latest.html) (Danecek et al., 2011), which is also well-documented and easy to use.


## 2. Download publicly available datasets used for the inference 
All of the annotations, genetic maps, polymorphism data, etc., are downloaded from publicly available datasets (see Appendix 1 Section 2 in [our paper](https://elifesciences.org/articles/76065) for details). Some of these, like the 100-vertebrate alignment MAF files used to generate phastCons scores and the 1000 Genomes Project VCF files, are quite large, so be sure you have ~ 1TB available if you plan to download everything at once. Many of the larger files (like the 100-vertebrate alignment MAF files) will only be used once to create smaller processed datasets, so the large originals can be deleted. Alternatively, some of the processing can be done by streaming remote files into a pipe using `curl`, although you will need to a stable connection for this option (see, for example, `/1000_genomes/sh/recodesubpop.sh` script, which processes a stream of 1000 Genomes Project VCF file data directly).

Most file paths used for input/output processing data are found in the `FileStruct` data structure defined in `lsm_python/classes/filestruct.py` and will be automatically updated once you set the `root_dir` variable (see above). Some file paths are also embedded in the relevant bash scripts or python files directly, so these must be updated accordingly (a search for `/Users/davidmurphy/`, `/ifs/data` and `ifs/scratch` should help you track down any leftover hard coded paths to update for in system, or better yet incorporate into the `FileStruct` data structure).

## 3. Pre-calculating background selection maps using `calc_bkgd`
The most intensive step in the pre-processing of data is the calculation of background selection maps using the modified version of Graham McVicker's `calc_bkgd` program (McVicker et al., 2009) found in `debug_bkgd/`. To speed up this step, we made a version of the code that calculates effects of background selection on partitions of chromosomes (we used 10M bp partitions in our results) that are created in parallel. The partitions are then reassembled into whole chromosome maps when they are all finishing running. Once the `calc_bkgd` exectutable file is compiled, you won't have to directly interact with the c program at all - the python program `lsm_python/precalc/calc_bkgd.py` will read a few parameters, built a configuration file and then use this to run the c program. 

The `lsm_python/precalc/calc_bkgd.py` program can be called by the shell script `sh/calc_bkgd.sh` by setting arguments in the command line. These arguments are: 
- `init`: `RunStruct` initialization file (see Quickstart Guide 3. Create new `RunStruct` file above)
- `chrom`: chromosome (chr1-chr22)
- `coef`: the selection coefficient for the effects of background selection given as a float (i.e., 10^-2 is input 0.01)
- `anno`: annotation label for annotation files of selected sites that will be used for calculating effects of background selection. Annotations files should be formatted and stored as they are currently found in the folder `/data/bsanno`. For example, the full path to annotations for chromosome 1 for the annotation `fish_cons94_new` would be: `<ROOT_DIR>/data/bsanno/fish_cons94_new/chr1.fish_cons94_new.bed`. `anno` only needs the label `fish_cons94_new` and the code will fill in the full path. If you introduce new background selection targets, make sure they conform to this formatting.
- `udel`: deleterious mutation rate (fixed at a default of 7.4e-8)
- `fexe`: path to executable c program `calc_bkgd`
- `pidx`: the partition index of the chromosome to calculate effects of background selection (units of 0:length(chrom) / partition size). We used a default partition size of 10M bp; for instance chr1, which is just under 250M bp, would have 25 10M bp partitions and `pidx` from 0 to 24
- `plen`: this is the above referenced partition size for breaking up chromosomes. We used 10M bp.
- `pts`: this is a flag (set to 1 for true, 0 for false) to indicate whether or not to apply a technical fix in the background selection calculation to avoid a rare loss of precision that occurs when the `calc_bkgd` program very rarely "jumps" over selected sites and misses their effect on neutral diversity in the resulting maps. In practice, setting this to 0 (false) will have very little effect on precision in most cases and will make the program run much faster. We set it to 1 (true) to deal with some of the issues of precision described in Appendix 1 Section 1.5 of the paper.

The key data files needed for calculating the effects of background selection are the set of annotations (set in the command line as describd above) and the genetic map. The genetic map is set in the `RunStruct` and the program `lsm_python/precalc/calc_bkgd.py` will look up the file paths to the genetic map files automatically. The files we used can be found in `data/genetic_maps/AA_Map_10kb` in the Google Drive folder linked above. 

NOTE: An additional wrapper script `sh/bkgd_batch.sh` calls `sh/calc_bkgd.sh` for a given range of parameters (for example, chromomosomes 1-10, all partition indices, a given range of selection coefficients, etc.). This can speed up the process of running dozens or hundreds of jobs in parallel instead of writing all of the commands to `sh/calc_bkgd.sh` independently.

The program `lsm_python/precalc/process_partitions.py` is used reassemble the partitioned chromosome maps once they are finished running. Once again, a shell wrapper `sh/prcs_bparts.sh` executes this code with arguments from the command line. These arguments are: 
- `chrom`: the chromosome
- `gmap`: name of the genetic map used (set in the `RunStruct` initilization file described above)
- `bdir`: name of the folder where background selection maps are stored (set in the `RunStruct` initilization file described above)
- `anno`: annotation label (as described in previous step)
- `coef`: selection coefficient (as described in previous step)
- `plen`: partition size (as described in previous step)

After running this script, a new background selection map made up of the merged partitions will be created and saved in the same directory as the partitioned files (`<ROOT_DIR>/precalc/<bdir>` where `<bdir>` is the argument referenced above and stored in the `RunStruct`; typically `bdir` and `anno` will be the same). These merged files are the raw input used in the data compression step described above in the quickstart guide prior to inference (4. Compress data).

## 4. Pre-calculating selective sweep maps
The calculation of selective sweep maps requires far less computation power and is therefore much more straightforward. The program `lsm_python/precalc/sweep_map.py` calculates these maps and is called by the shell wrapper `sh/sweep_map.sh`. The arguments given in the command line are:
- `chrom`: the chromosome (chr1-chr22)
- `anno`: the annotation label. As with the `anno` argument used in calculating background selection maps, the annotation label will be used to fill in the full file paths by the code. For example, the full path to nonsynonymous substitutions on chromosome 1 would be: `<ROOT_DIR>/data/csanno/nonsyn/chr1.nonsyn.npy`. `anno` only needs the label `nonsyn` and the code will fill in the rest. If you introduce new selective sweep targets, make sure they conform to this formatting.
- `coef`: the selection coefficient as a float (see `coef` for background selection maps above).
- `lim`: an optional flag to ignore substitutions beyond a certain limit where they are estimated to have no effect on diversity levels for a given neutral site. Setting this to "True" can speed up the calculations somewhat but in practice it doesn't make a very big difference because the calculation of sweep maps is very efficient. 

As with the maps of background selection described above, the key data files needed for calculating the effects of selective sweeps are the set of substitutions (set in the command line as describd above) and the genetic map. The genetic map files can be found in `data/genetic_maps/AA_Map_10kb` in the Google Drive folder linked above.

Completed maps will be saved in the path (`<ROOT_DIR>/precalc/<cdir>` where `<cdir>` is set in the `RunStruct`; typically `cdir` and `anno` will be the same).

### NOTE: This documentation is being expanded to include additional use cases. If you have questions beyond what is presented here, please email David Murphy (murphyd@omrf.org).

