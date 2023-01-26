#!/bin/bash -l

# Runs synthesis pipeline
# Inputs are lambda, optimal angle, discretization size, and index to design coupling into
# Specific for 45 CLO
#
# Example
#	qsub run_synth_bash_45CLO.sh 1300 10 shallow up 5
#		The above runs the synthesis for wavelength of 1300nm, angle 10 deg, partial etch, upwards coupling, 5 um disc
#		for coupling into air

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
matlab -nodisplay -r "addpath('/project/siphot/bz/code/utility'); addpath('/projectnb/siphot/bz/git/grating_designs/45CLO_2022_12_X'); f_run_designspace_duallevel_45CLO( $1, $2, '$3', '$4', $5 ); exit"