#!/bin/bash -l

# Runs synthesis pipeline
# Inputs are lambda, optimal angle, discretization size, 'up' or 'down', geometry ('a', 'b', etc)
# Specific for aim grating
#
# Example
#	qsub run_synth_bash_aim.sh 1300 10 5 up 'a'
#		The above runs the synthesis for wavelength of 1300nm, angle 10 deg, and disc. size of 10nm
#		for optimizing unidirectionality upwards
#       design type 'a' (default)

# set job name
#$ -N grating_synth

# set # cores to use
#$ -pe omp 28

# specify hard time limit for job
#$ -l h_rt=23:00:00

# merge output and error files in one
#$ -j y

# send me email if job finishes or aborted
#$ -m ea

# load matlab
module load matlab/2019a

# run script
# matlab -nodisplay -r "addpath('/project/siphot/bz/code/utility'); addpath('/projectnb/siphot/bz/git/grating_synthesis/chips/210716_AIM_SOPAs'); f_run_designspace_2level_aim( $1, $2, $3, '$4', '$5' ); exit"
matlab -nodisplay -r "addpath('/projectnb/siphot/howard/git/grating_designs/aim_2021_07_16_SOPAs'); f_run_designspace_2level_aim( $1, $2, $3, '$4', '$5' ); exit"