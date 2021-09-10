This README file explains the basics of how to set up and run this code on a new local machine. Everything was run using Python 2.7, so some compatibility changes may be needed for Python 3.x+. The basic steps are: 1. downloading and installing all software, 2. updating file paths in python and bash scripts, 3. downloading public datasets used in analysis, 4. formatting and processing data, 5. running the inference, 6. analyzing the inference results

-

# 0. Using B-maps

If you are just interested in using the B-maps that are referenced in our paper, these can be found in the `/Bmaps/` directory, along with relevant README files. These are formatted using the hg19 genome build.

-

# 1. Downloading and installing all software
First, download the entire contents of this repository, including all python code, shell scripts and the modified version of McVicker's calc_bkgd program (McVicker et al., 2009) in `/debug_bkgd/`. You will probably need to install some dependencies to run calc_bkgd, depending on your system. This step can bit a bit of a pain, depending on your programming skill level. I have included a list of the required C libraries and some notes on installation that were helpful in my case in the `notes-on-dependencies.txt` file in the `/debug_bkgd/` folder.

Additional software that are used for data processing are `phastCons` and `phyloFit` from the PHAST suite of tools (Siepel and Haussler 2004; Siepel, Bejerano et al. 2005). 

-

# 2. Updating file paths in python and bash scripts

Both python and bash scripts have a number of embedded paths that are used to call either other scripts or for loading data. The vast majority of these in python are linked to the `root_dir` variable in the `/classes/runstruct.py` file, so updating this path to fit your own directory structure will carry through 99% of downstream paths. There is a small subset of paths to `/ifs/scratch/c2b2/gs_lab/dam2214/scratch` that you will need to search+replace to match your directory structure. 

Similarly, the bash scripts stored in the `/sh/` folder generally include paths to the python binary file from my system `py=/path/` and calls to the python program that will be called by the shell script `prog=/path/` that will need to be updated to match the current file structure.



