#!/bin/bash
#SBATCH --partition       h3c
#SBATCH --time            12:00:00
#SBATCH --nodes           1
#SBATCH --ntasks-per-node 24
#SBATCH --cpus-per-task   1
#SBATCH --error           out_%j.err
#SBATCH --output          out_%j.log

################################################################################
ulimit -s unlimited

# module load vasp/5.4.4
# module load vasp/5.4.4_soc
module load vasp/6.1.0
# module load vasp/6.1.0_omp

################################################################################

#export I_MPI_PMI_LIBRARY=libpmi2.so

echo "============================================================"
module list
env | grep "MKLROOT="
echo "============================================================"
echo "Job ID: $SLURM_JOB_NAME"
echo "Job name: $SLURM_JOB_NAME"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "Number of processors: $SLURM_NTASKS"
echo "Task is running on the following nodes:"
echo $SLURM_JOB_NODELIST
echo "OMP_NUM_THREADS = $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo

################################################################################
# Modified For LOOP
################################################################################
START=1
END=2000
NDIGIT=$(echo -n ${END} | wc -c)
RUNDIR=${PWD}
PREFIX="run"

for i in `seq ${START} ${END}`
do
    (( j = i - 8 ))
    ii=`printf "%0${NDIGIT}d" ${i}`
    jj=`printf "%0${NDIGIT}d" ${j}`

    if [[ -d "${PREFIX}/${ii}" ]]; then

        # goto the working directory
        cd ${PREFIX}/${ii}

        # Skip this directory if job is running or has ended
        if [[ -f RUNNING || -f ENDED ]]; then
            # Do nothing, go back to the parent directory
            cd ${RUNDIR}
            continue
        fi

        # Mark the directory as RUNNING
        touch RUNNING
        echo "#### RUNNING in DIR: ${RUNDIR}/${PREFIX}/${ii}"
        sleep 1

        # Copy the CHGCAR from the nearest finished job to accelerate current job
        [[ -s ../${jj}/CHGCAR && -f ../${jj}/ENDED ]] && cp ../${jj}/CHGCAR .

        ############################################################
        # Run VASP
        ############################################################
        # srun vasp_std
        srun vasp_gam
        # srun vasp_ncl
        ############################################################

        # Check if jobs ended correctly
        if grep 'Total CPU' OUTCAR >& /dev/null; then
            # If so, mark the job as ENDED
            touch ENDED
        else
            rm ENDED 2> /dev/null
        fi

        # Delete redundant files
        rm RUNNING CHG vasprun.xml >& /dev/null

        # Go back to the parent directory
        cd ${RUNDIR}
    fi
done
