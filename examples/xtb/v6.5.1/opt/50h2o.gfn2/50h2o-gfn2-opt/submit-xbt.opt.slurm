#!/bin/bash
#SBATCH --job-name=50h2o-gfn2-opt
#SBATCH --output=50h2o-gfn2-opt.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1-00:00:00
#SBATCH --cluster=smp

cd $SLURM_SUBMIT_DIR
module purge

# Sets the path to your xtb excutible
xtb=/ihome/jkeith/amm503/codes/xtb-6.5.1/bin/xtb

# Recommendations from https://xtb-docs.readthedocs.io/en/latest/setup.html#setting-up-xtb
ulimit -s unlimited
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# System properties
xyz_file='./50h2o.pm.xyz'
charge=0  # Charge of the system.
multiplicity=1  # Specifies multiplicity of the system.
spin_num=$multiplicity-1

# Calculation properties
gfn_level=2  # 0, 1, or 2. 2 is generally prefered. For more information: https://xtb-docs.readthedocs.io/en/latest/commandline.html#options
runtype=opt
conv_level=normal  # loose, lax, normal, tight, vtight, or extreme optimization convergence. For more information: https://xtb-docs.readthedocs.io/en/latest/optimization.html#optimization-levels
opt_cycles=500  # Maximum number of optimization cycles.





# Moves calculation to SLURM scratch.
cp $SLURM_SUBMIT_DIR/$xyz_file $SLURM_SCRATCH

# Sets a trap to copy over files when job ends or is canceled.
run_on_exit(){
    cp -r $SLURM_SCRATCH/xtbopt.log $SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.log
    cp $SLURM_SCRATCH/xtbopt.xyz $SLURM_SUBMIT_DIR/$SLURM_JOB_NAME.xyz
}
trap run_on_exit EXIT 

# Runs xtb optimizations.
cd $SLURM_SCRATCH
if [ $multiplicity -gt 1 ]
then
    # Uses unrestricted Hartree Fock if there are any unpaired electrons.
    $xtb $xyz_file --$runtype $conv_level --gfn $gfn_level --charge $charge --uhf $spin_num --cycles $opt_cycles
else
    # Uses restricted Hartree Fock for singlet states. 
    $xtb $xyz_file --$runtype $conv_level --gfn $gfn_level --charge $charge --cycles $opt_cycles
fi
