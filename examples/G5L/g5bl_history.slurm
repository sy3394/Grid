#!/bin/bash -l
#SBATCH --job-name=g5bl_hist
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=04:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -o g5bl_hist_%j.out

# ---------------------------------------------------------------------------
# Performance analysis on 32^4: raw Euclidean residual & #converged-modes vs Krylov
# dimension for the two Grid solvers on the SAME shift-inverted operator
# (D_W - sigma)^{-1}:
#   refinedarnoldi - RefinedArnoldi (single-vector Arnoldi + refined Ritz)
#   g5bl           - Gamma5BlockLanczos (gamma5-metric block method)
# No SLEPc needed.
#
# Writes <out>.{refinedarnoldi,g5bl}.hist with columns
#   krylov_dim  min_raw_res  n_below_1e-6  n_below_1e-10
# Krylov dim = number of (D_W - sigma)^{-1} applications.
#
# Pull the .hist files back to ~/BNL/G5BL/data and plot with g5bl_performance.ipynb.
# The --compare run prints a 3-row head-to-head (g5bl / arnoldi / blockArnoldi).
# ---------------------------------------------------------------------------

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 7 6 5 4)
export NUMA_MAP=(3 3 1 1 2 2 0 0)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUM=\${NUMA_MAP[\$SLURM_LOCALID]}
export HIP_VISIBLE_DEVICES=\$GPU
unset ROCR_VISIBLE_DEVICES
exec numactl -m \$NUM -N \$NUM \$*
EOF
chmod +x ./select_gpu

root=$HOME/Frontier/Grid/systems/Frontier
source ${root}/sourceme.sh
export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1

EXE=./Example_gamma5_block_lanczos
CONFIG=$HOME/path/to/snapshots/702/U_smr_3.000000   # <-- EDIT
OUT=hist702                                          # tag for the .hist files

SIGMA=0.78        # centre of the reference cluster Re in [0.749,0.812]
MASS=0
STEPS=40          # RefinedArnoldi Krylov dim = 2*STEPS = 80; g5bl block steps = 40
KRYLOV=80
STOL=1e-11        # inner CG tolerance (double precision)
SITER=30000

GRIDOPT="--accelerator-threads 8 --comms-overlap --shm 4096 --shm-mpi 0"

# Residual-vs-Krylov-dim history for RefinedArnoldi and Gamma5BlockLanczos.
srun ./select_gpu ${EXE} --grid 32.32.32.32 --mpi 1.2.2.2 ${GRIDOPT} \
    --config ${CONFIG} --mass ${MASS} --shift ${SIGMA} \
    --stol ${STOL} --siter ${SITER} --steps ${STEPS} --krylov ${KRYLOV} \
    --history --out ${OUT}

# (optional) one-shot head-to-head summary at equal Krylov dim:
srun ./select_gpu ${EXE} --grid 32.32.32.32 --mpi 1.2.2.2 ${GRIDOPT} \
    --config ${CONFIG} --mass ${MASS} --shift ${SIGMA} \
    --stol ${STOL} --siter ${SITER} --steps ${STEPS} --krylov ${KRYLOV} --accept 1e-6 --compare

echo "Done.  Copy ${OUT}.{refinedarnoldi,g5bl}.hist to ~/BNL/G5BL/data and run g5bl_performance.ipynb"
