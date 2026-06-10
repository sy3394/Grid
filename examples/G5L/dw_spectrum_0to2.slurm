#!/bin/bash -l
#SBATCH --job-name=DWspec_702
#SBATCH --partition=extended
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=10:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --out=out/DWspec_702-%j.out

# ---------------------------------------------------------------------------
# Wilson D_W eigenvalues in the real-axis window Re(lambda) in [LO,HI], via a
# shift-invert sweep (Wilson_DW_spectrum; RefinedArnoldi by default, --g5bl for
# gamma5-block Lanczos).  One config per job; output rows  <tag>  Re(lambda)  Im(lambda)
# (sorted by |Im|, near-real first).  Resubmit per config; narrow WINDOW for the cluster.
#
# 4 nodes x 8 GCDs = 32 ranks; 32^4 with --mpi 4.2.2.2.
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
EXE=${root}/examples/Wilson_DW_spectrum
CONFIG=${root}/HMC/32cube-rho0.124-tau4/ckpoint_EODWF_lat_smr.702
OUT=evals_DW_702.dat
TAG=3.000000                                          # MD-time label (col 1)

# ---- window + solver parameters ----
MASS=0                # reference used m = 0
WINDOW=0:2            # Re(lambda) window (narrow to the cluster with e.g. 0.70:0.85)
NSHIFT=11             # real shifts tiling the window (validation; raise to 21 for the full run)
KRYLOV=80             # RefinedArnoldi Krylov dimension per shift
ACCEPT=1e-8           # keep modes with raw residual < ACCEPT
DEDUPE=1e-6           # merge modes within this distance in C
STOL=1e-9             # inner CG tolerance: ~10x tighter than ACCEPT is plenty.  Shifts
                      # that land on the cluster make (D_W-sigma) near-singular, so 1e-11
                      # would cost ~8500 iters/solve for no gain in the final eigenvalues.
SITER=30000           # inner CG max iterations
# g5bl-only knobs (used when --g5bl is appended below):
STEPS=40              # gamma5-block Lanczos block steps (Krylov dim = 2*STEPS)
WANTED=12             # wanted conjugate pairs per window
CYCLES=6              # thick-restart cycles per shift


vol=32.32.32.32
mpi=4.2.2.2
PARAMS="--mpi $mpi --accelerator-threads 32 --shm 2048 --shm-mpi 0 --grid $vol --comms-overlap" #--comms-sequential

srun ./select_gpu ${EXE} $PARAMS \
    --config ${CONFIG} --mass ${MASS} \
    --window ${WINDOW} --nshift ${NSHIFT} --krylov ${KRYLOV} \
    --accept ${ACCEPT} --dedupe ${DEDUPE} --stol ${STOL} --siter ${SITER} \
    --tag ${TAG} --out ${OUT}
# To produce the spectrum with gamma5-block Lanczos instead, append:
#     --g5bl --steps ${STEPS} --wanted ${WANTED} --cycles ${CYCLES}

# Output: ${OUT} = D_W eigenvalues with Re(lambda) in the window, de-duplicated.
# Compare:  python3 G5L/compare_evals.py evals_DW_702_2.dat ${OUT} --mdtime 3.0
