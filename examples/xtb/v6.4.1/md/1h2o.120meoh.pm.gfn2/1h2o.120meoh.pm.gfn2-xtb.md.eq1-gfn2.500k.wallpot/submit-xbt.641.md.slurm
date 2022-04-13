#!/bin/bash
#SBATCH --job-name=1h2o.120meoh.pm.gfn2-xtb.md.eq1-gfn2.500k.wallpot
#SBATCH --output=1h2o.120meoh.pm.gfn2-xtb.md.eq1-gfn2.500k.wallpot.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=3-00:00:00
#SBATCH --cluster=smp

##### Specify for each job #####

# System properties
# Should already be optimized with xtb.
xyz_file='1h2o.120meoh.pm.gfn2.xyz'

charge=0  # Charge of the system.
multiplicity=1  # Specifies multiplicity of the system.

# Writes xtb MD input file.
# For more information:
# https://xtb-docs.readthedocs.io/en/latest/md.html#parameters
## Note: Specifying a temperature sets nvt to true. No temperature makes it NVE.
input_file='1h2o.120meoh.pm.gfn2-xtb.md.eq1-gfn2.500k.wallpot.inp'  # Must end in '.inp'.
cat > $input_file << 'EOF'
$md
    restart  = false
    time     =  10.0  # in ps (1000 fs = 1 ps; 1 fs = 0.001 ps)
    step     =   1.0  # in fs
    dump     =  50.0  # in fs
    temp     = 500.0  # in K
    nvt      = true
    velo     = false
    hmass    =   1
    shake    =   0
    sccacc   =   1.0
$end

$wall
   potential = logfermi
   sphere: 23.62158, all
   sphere: 2.83459, 1,2,3
   temp = 300.0
   beta = 6
$end

$write
   wiberg=false
   dipole=false
   charges=false
   mulliken=false
   orbital energies=false
   inertia=false
   distances=false
   final struct=false
   geosum=false
$end
EOF

# Calculation properties
runtype=md
gfn_level=2  # 0, 1, or 2. 2 is generally prefered.

##### Starts job #####

# Start every calculation with no preloaded packages.
cd $SLURM_SUBMIT_DIR
module purge

# Sets the path to your xtb excutible
xtb=/ihome/jkeith/amm503/codes/xtb-6.4.1/bin/xtb

# Recommendations from
# https://xtb-docs.readthedocs.io/en/latest/setup.html#setting-up-xtb
ulimit -s unlimited
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Moves everything to scratch directory to run.
cp $SLURM_SUBMIT_DIR/$input_file $SLURM_SCRATCH
cp $SLURM_SUBMIT_DIR/$xyz_file $SLURM_SCRATCH
cp $SLURM_SUBMIT_DIR/mdrestart $SLURM_SCRATCH 2>/dev/null || :

# Sets a trap to copy over files when job ends or is canceled.
run_on_exit(){
    cp $SLURM_SCRATCH/* $SLURM_SUBMIT_DIR/
}
trap run_on_exit EXIT 

# Runs xtb optimizations.
cd $SLURM_SCRATCH
if [ $multiplicity -gt 1 ]
then
    # Uses unrestricted Hartree Fock if there are any unpaired electrons.
    # Unsure if this is right (not checked).
    spin_num=$multiplicity-1
    $xtb $xyz_file --$runtype --input $input_file --gfn $gfn_level --charge $charge --uhf $spin_num --verbose
else
    # Uses restricted Hartree Fock for singlet states. 
    $xtb $xyz_file --$runtype --input $input_file --gfn $gfn_level --charge $charge --verbose
fi
