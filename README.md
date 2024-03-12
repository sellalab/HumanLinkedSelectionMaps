# OVERVIEW
This README file explains the basics of how to set up and run the `HumanLinkedSelectionMaps` code and replicate the results presented in [our paper](https://elifesciences.org/articles/76065). The code can be adapted to work with data from different species or different human datasets not explored in our paper (see "Advanced usage" section). If you are just interested in downloading and using the B-maps, these can be found in the `Bmaps/` directory, along with relevant README files. B-maps are formatted using the hg19 genome build.

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
First, download the entire contents of this repository. Name the folder `lsm_python` (or whatever you like). Paths to input and output files are set in the `RunStruct` class, which indlues a hard-coded `root_dir` variable at the top of the file `classes/runstruct.py`. This was set to have two options, a local option and a remote option, depending on the User path on your system. Both paths must be updated here:

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
The `RunStruct` class includes a representation of the complete file structure used for input and output as well as inference settings, inference results, and other important variables. The code for this class is found in `classes/runstruct.py` (see above section on the root directory). For a given configuration of the inference, a new `RunStruct` is initialized and saved to a text file that is used as the input to a number of other programs. The file `likelihood/initialize_runstruct.py` is set up to initialize a `RunStruct` and save it in text format for the inference results presented as "CADD 6%" in the main paper. It will also create some folders for initialization and results files. Run this code and check that the file `/<ROOT_DIR>/result/init_files/YRI.cadd94_gmask_v1.6_without_bstat.BS1.6.CS0.0.NOT_STARTED.initial.txt` was created in your root directory. 

To initialize a `RunStruct` for the other main result shown in the paper (phastCons 6%), change the line `b_an_1 = 'cadd94_gmask_v1.6_without_bstat'` to `b_an_1 = 'fish_cons94_new'` in `likelihood/initialize_runstruct.py` and run the code again. This should create a comparable text file with the new annotation in the file name.

## 4. Compress data
As detailed in Appendix 1 Section 1.3, we compress neutral sites into "bins" where the values of exogenous parameters or pre-calculated "raw" B values are constant, which vastly downscales the number of calculations needed in the inference without any loss of precision in our predictions. The file `precalc/lh_inputs.py` will perform this step by taking as input the `RunStruct` file created in the previous step. The processing pipeline takes one chromosome at a time and can be used with external command line arguments to run each autosome in parallel on a cluster (see code segments in the `main` function that are commented out). In the current configuration of the `main` function, a `for` loop will just go through each chromosome in serial. Note that this step may take 10-20 minutes (depending on your system).

## 5. Run inference
As described in Appendix 1 Section 1.4, we run the inference in two optimization steps. In the first step, we run the optimization from 15 different starting conditions and keep the 3 results with the greatest likelihood. In the second step, we take the mean parameters of these 3 "best" runs from stage 1 and use these as the input parameters to the optimization. 

To run the first step, run the script `likelihood/inf_step_1.py`. This will read the `RunStruct` file (created above) and loop through the 15 different starting conditions (indexed 0-14), leading to 15 different results files in the folder `<ROOT_DIR>/result/final_files/cadd94_gmask_v1.6_without_bstat/`, with index 0-14 in the file labels. Note that this step is much faster if you can run each of the 15 optimizations in parallel and take the index as an external command line argument (see code segments in the `main` function that are commented out).

After all 15 optimizations finish, step 1 is complete and you can move on to step 2. To do so, run the script `likelihood/inf_step_2.py`. As in the above steps, this script is hard coded to continue running through the example for CADD 6% which we've started in the previous steps. The `main` function will just read the name of the results folder (i.e., `cadd94_gmask_v1.6_without_bstat`) as input. After the step 2 optimization is complete, a new file with the suffix `.composite.txt` will be saved to the folder `<ROOT_DIR>/result/final_files/cadd94_gmask_v1.6_without_bstat/`. This file contains the best-fitting parameters the current model configuration and data.

## 6. Create figures
After the inference is run and the final `.composite.txt` file has been created, you can use the inferred parameters to perform the analyses presented in the paper. Here we will cover the figures from the main paper (Figs. 2-5); to generate figures from the Appendix, see **ADVANCED USAGE** below.

Each of the figures requires an analysis step that results in a small table of values that are used to render the plots in Figs. 2-5. To simplify, each of these steps has been put into a single python file: `figures/main_figs_analysis.py`. Running this script will run each analysis for Figs. 2-5 and save the results in text or numpy array format to the folder `cadd94_gmask_v1.6_without_bstat`. You can explore each of the called analysis functions in more detail if you'd like to better understand the analyses, or change certain paramters (e.g., you could change the chromosome in figure 2a from chr1 to chr5).

To plot the figures, use the script `figures/main_figs_new_computer.py`. This script is set up like a python notebook with individual cells for each figure. The individual cells will make each of the plots (Figs. 2-5). You can also just run the whole script and generate all of the figures. Saved files will be written to `<ROOT_DIR>/result/final_files/mainfigs`.

**NOTE: The `results/final_files/` folder includes some additional results folders that are used for mutation rate estimates plotted in figure 4.**

### NOTE: Code documentation is being updated to include a detailed guide for more advanced usage. If you have any immediate questions, please email David Murphy (murphyd@omrf.org).
