log
-------------
# (02/04/2025): fixed a bug that was removing much of the kinematic data across 30 subjects and was potentially corrupting other subject kinematic data. all data after this date should be goood according to validation figures.

# (02/06/2025): fixed another bug in the PSD generation (fixed baselines designated with "new") where some of the baselines were not being layered correctly and therefore multiple baselining procedures did not work. 

# (06/17/2025): Some running notes from .sh scripts

## to assign certain libraries for specific R instances, but this doesn't seem to work?
mkdir ~/R/x86_64-pc-linux-gnu-library/4.4
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.4

## loading R... Just make the conda enviornment
module load R/4.4
rm -f .Rprofile

(05/08/2025) JS, this fixes an error where Rmpi isn't able to load upon starting R
this will need to most likely be reinitiated. This can be done using VSCode or MobaXTerm (and others?).

other loading options for rmpi runs:
module load R/4.4
module load gcc/12.2.0 openmpi/4.1.6 rmpi/4.4
ml gcc/12.2.0 openmpi/4.1.6 R/4.3
ln -s /home/jsalminen/ .Rprofile
ln -s /apps/rmpi/conf/Rprofile .Rprofile

(05/09/2025) JS, kept getting an error where the R instance couldn't detect rmpi, so I had to contact IT and they remade the rmpi instance for me.

## MAKE R CONDA ENVIORNMENT
(06/17/2025) JS, after consulting with IT many times this is a guaranteed solution:
you can build your own R using conda-forge ...

ml conda
conda create -n rmpi-env r-rmpi
conda activate rmpi-env

NOTE: use "conda deactivate" to deactivate the environment in the terminal

### install packages after loading R.
this takes a very long time sadly:

install.packages('lifecycle','igraph','genlasso','gtable','scales','pillar','lmerTest','emmeans','R.devices','dplyr','tibble','Matrix','lattice','ggplot2','gridExtra')
install.packages('backports','pbkrtest');

NOTE: make sure to get all co-dependencies. Hipergator doesn't automatically do this?

### to use env
add path to each batch script:
export PATH=/blue/mygroup/$USER/conda/envs/project1/bin:$PATH (e.g., "PATH=/blue/dferris/$USER/.cond/envs/rmpi-env/bin:$PATH")

OR

ml conda
conda activate name_of_env (e.g., conda activate rmpi-env)

## GENERIC PRINTS

HPC_R_DIR - installation directory
HPC_R_BIN - executable directory
HPC_R_LIB - library directory
HPC_R_INCLUDE - includes directory

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "library directory: $R_LIBS_USER"
echo "installation directory: $HPC_R_DIR"
echo "executable directory: $HPC_R_BIN"
echo "library directory: $HPC_R_LIB"
echo "includes directory: $HPC_R_INCLUDE"