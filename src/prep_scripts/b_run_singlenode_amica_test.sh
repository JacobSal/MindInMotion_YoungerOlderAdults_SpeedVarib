#!/bin/bash
#SBATCH --job-name=ALL_AMICA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsalminen@ufl.edu
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_slurm_logs/%j_amica_out.log
#SBATCH --nodes=64
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000mb
#SBATCH --distribution=cyclic:cyclic
#SBATCH --time=02:00:00
#SBATCH --account=dferris
#SBATCH --constraint=el9        # el9 is REHL 9 login, el8 has the original coding on it
#SBATCH --qos=dferris-b
#SBATCH --partition=hpg-default
# NOTE: (04/22/2023) SalminenJ, Seems to time out after ~15 subaject runs (6hr time limit). Moving to a 48hr cycle ( I think this is max for hpg2-compute).
# sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/b_run_singlenode_amica_test.sh
# (07/17/2025) JS, updated to use newer compiler and libraries. This is a test script for the new amica17ub binary.
# had to make some formatting edits to amica17.f90: updated random_seed use to force single integer & updated compile script to
# modernize to HPG rehl9 node configuration. I removed the using of intel/2020 libraries in favor of gcc/12.2.0 and openmpi/5.0.7
# going to pull request EEGLAB's amica lib to see if I can't get an update passed.

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "slurm_mem_per_cpu $SLURM_MEM_PER_CPU"
echo "slurm_mem_per_gpu $SLURM_MEM_PER_GPU"
echo "slurm_mem_per_node $SLURM_MEM_PER_NODE"

module purge
module load gcc/12.2.0 openmpi/5.0.7 openblas lapack mkl

#%% PARAMS
export DONOT_RECREATE=false;

export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej"

export SUBJ_RUN=("H1002" "H1004")

#%% LOOP
for s in ${SUBJ_RUN[@]}
do
	export param_f=$SUBJ_DIR/$s/clean/*.param
	export chk_f=$SUBJ_DIR/$s/clean/W
	echo "Processing $s eeg file..."
	if test -f "$chk_f" && $DONOT_RECREATE;
	then
		echo "ICA weight file is already generated: $chk_f"
	else
		if test -f $param_f;
		then
			echo "Calculating ICA weights..."
			srun --mpi=pmix_v5 /blue/dferris/share/s.peterson/rhel9_amica_15/amica17ub $param_f			
			#%% This "wait" may be needed to ensure shared libraries aren't accessed multiple times?
			wait
			echo "done: $s"
		else
			echo "$param_f does not exist. Preprocess EEG for $s."
		fi
	fi
done

