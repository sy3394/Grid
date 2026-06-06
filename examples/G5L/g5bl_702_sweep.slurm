#!/bin/bash -l
#SBATCH --job-name=g5bl_702
#SBATCH --partition=extended
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
##SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=8
#SBATCH --time=24:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --out=out/g5bl_702-%j.out

# ---------------------------------------------------------------------------
# gamma5-Block Lanczos with a SHIFT-INVERT SWEEP over real sigma, in ONE job.
# Resolves the near-real interior D_W spectrum of config 702 band by band:
# each sigma converges the eigenvalues nearest it (lambda = sigma + 1/theta),
# and the code accumulates + de-duplicates across windows.
#
# Reference cluster (evals_DW_702_2.dat, MD-time slice): Re in [0.749, 0.812].
# The default sweep 0.70:0.85:8 tiles that band; widen SWEEP for more of the
# interior real axis (e.g. 0.0:1.4:15).
#
# 4 nodes x 8 GCDs = 32 ranks; 32^4 with --mpi 2.2.2.4.
# ---------------------------------------------------------------------------

cat << EOF > ./select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 4 5 6 7)
export NUMA_MAP=(3 3 1 1 0 0 2 2)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUM=\${NUMA_MAP[\$SLURM_LOCALID]}
export ROCR_VISIBLE_DEVICES=\$GPU
echo RANK \$SLURM_LOCALID using GPU \$GPU NUMA \${NUM}
exec numactl -m \$NUM -N \$NUM \$*
EOF
chmod +x ./select_gpu

mkdir -p out

# Grid Frontier environment (adjust to your tree)
root=$HOME/Frontier/Grid/systems/Frontier
source ${root}/sourceme.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1

# ---- paths (EDIT) ----
EXE=./Example_gamma5_block_lanczos
CONFIG=$HOME/path/to/snapshots/702/U_smr_3.000000     # <-- your snapshot
OUT=g5bl_702_md3.dat
TAG=3.000000                                          # MD-time label (col 1)

# ---- sweep + solver parameters ----
SWEEP=0.70:0.85:8     # n real shifts tiling the reference cluster
MASS=0                # reference used m = 0
STEPS=40              # Lanczos steps per cycle (active Krylov 2*STEPS)
WANTED=12             # wanted conjugate pairs per shift window
CYCLES=6              # thick-restart cycles per shift
TOL=1e-7              # Lanczos (theta-space) tolerance
ACCEPT=1e-4           # keep modes with theta-residual < ACCEPT
DEDUPE=1e-5           # merge modes within this distance in C
STOL=1e-8             # inner normal-equations CG tolerance
SITER=20000           # inner CG max iterations

# NUMA/CPU binding is handled by numactl inside select_gpu; do NOT pass
# --cpu-bind=map_ldom to srun (Frontier rejects it unless whole-node, and it
# conflicts with the numactl binding).
srun ./select_gpu ${EXE} \
    --grid 32.32.32.32 --mpi 2.2.2.4 \
    --accelerator-threads 8 --comms-overlap --shm 4096 --shm-mpi 0 \
    --config ${CONFIG} --mass ${MASS} \
    --shift-sweep ${SWEEP} \
    --steps ${STEPS} --wanted ${WANTED} --cycles ${CYCLES} --tol ${TOL} \
    --accept ${ACCEPT} --dedupe ${DEDUPE} --stol ${STOL} --siter ${SITER} \
    --tag ${TAG} --out ${OUT}

# Output: ${OUT} = all converged, de-duplicated D_W eigenvalues across the sweep.
# Compare:  python3 G5L/compare_evals.py evals_DW_702_2.dat ${OUT} --mdtime 3.0
