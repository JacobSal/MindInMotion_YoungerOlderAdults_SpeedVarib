#!/bin/bash
#SBATCH --job-name=STAT_FLASSO # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks
#SBATCH --cpus-per-task=32 # Number of tasks on each node
#SBATCH --mem-per-cpu=8gb# Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=08:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_slurm_logs/%j_flasso_lmer_stat_nocb.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --constraint=el8 # el9 is REHL 9 login, el8 has the original coding on it
#SBATCH --partition=hpg-default # cluster to run on, use slurm command 'sinfo -s'
# sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/r_scripts/sbs_lme_mods/run_tf_flasso_lmer_stat_speed_nocb.sh

module purge
module load conda
conda activate rmpi-env

#--
# HPC_R_DIR - installation directory
# HPC_R_BIN - executable directory
# HPC_R_LIB - library directory
# HPC_R_INCLUDE - includes directory

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

# set linux workspace
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    TMP_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    TMP_PATH=$(realpath $0)
fi

export SCRIPT_DIR=$(dirname $TMP_PATH)
export STUDY_DIR=$SCRIPT_DIR
export SRC_DIR=$SCRIPT_DIR
cd $STUDY_DIR

# mkdir -p $SCRIPT_DIR/tmp
# export TMPDIR=$SCRIPT_DIR/tmp

# Create a temporary directory on scratch
mkdir -p $STUDY_DIR/_slurm_scratch/$SLURM_JOB_ID

# Kick off R
Rscript $SCRIPT_DIR/tf_flasso_lmer_stat_speed_nocb.R

# Cleanup local work directory
rm -rf $STUDY_DIR/_slurm_scratch/$SLURM_JOB_ID
