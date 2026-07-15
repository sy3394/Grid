#!/bin/bash -l
#SBATCH --job-name=solver_cmp
#SBATCH --partition=extended
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=06:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --out=out/compare_solvers-%j.out

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
# The --compare run prints a head-to-head at equal Krylov dim.
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

root=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier
source ${root}/sourceme.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1

# ---- paths (EDIT) ----
EXE=${root}/examples/Example_gamma5_block_lanczos
CONFIG=${root}/HMC/32cube-rho0.124-tau4/ckpoint_EODWF_lat_smr.702
OUT=hist702                                          # tag for the .hist files

SIGMA=0.78        # centre of the reference cluster Re in [0.749,0.812]
MASS=0
STEPS=40          # RefinedArnoldi Krylov dim = 2*STEPS = 80; g5bl block steps = 40
KRYLOV=80
# Loose inner CG: sigma=0.78 sits ON the cluster so (D_W-sigma) is near-singular and
# CG to 1e-11 costs ~8500 iters/solve.  Convergence is judged on the RAW D_W residual
# (and both solvers share the inner solve), so 1e-7 leaves the comparison curves intact
# while roughly halving the iterations.  The curves then plateau near ~1e-6.
STOL=1e-7
SITER=30000


vol=32.32.32.32
mpi=1.2.2.2       # 1 node x 8 GCDs
PARAMS="--mpi $mpi --accelerator-threads 32 --shm 2048 --shm-mpi 0 --grid $vol --comms-overlap" #--comms-sequential

# Residual-vs-Krylov-dim history for RefinedArnoldi and Gamma5BlockLanczos.
srun ./select_gpu ${EXE} $PARAMS \
    --config ${CONFIG} --mass ${MASS} --shift ${SIGMA} \
    --stol ${STOL} --siter ${SITER} --steps ${STEPS} --krylov ${KRYLOV} \
    --history --out ${OUT}

# (optional) one-shot head-to-head summary at equal Krylov dim:
srun ./select_gpu ${EXE} $PARAMS \
    --config ${CONFIG} --mass ${MASS} --shift ${SIGMA} \
    --stol ${STOL} --siter ${SITER} --steps ${STEPS} --krylov ${KRYLOV} --accept 1e-6 --compare

echo "Done.  Copy ${OUT}.{refinedarnoldi,g5bl}.hist to ~/BNL/G5BL/data and run g5bl_performance.ipynb"
