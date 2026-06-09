#!/bin/bash -l
#SBATCH --job-name=DWspec_0to2
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=12:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --out=DWspec_0to2-%j.out

# ---------------------------------------------------------------------------
# Wilson D_W eigenvalues in the real-axis window Re(lambda) in [0, 2], for the
# 32^4 rho0.124-tau4 ensemble, via a shift-invert sweep (RefinedArnoldi).
#
# SUBMIT FROM the configuration directory, e.g.
#   .../systems/Frontier/HMC/32cube-rho0.124-tau4/
#       sbatch <path-to-Grid>/examples/G5L/dw_spectrum_0to2.slurm
# so that $SLURM_SUBMIT_DIR holds the gauge configs; one output file per config:
#   evals_DW_<id>.dat   with rows   <id>  Re(lambda)  Im(lambda).
#
# 4 nodes x 8 GCDs = 32 ranks; 32^4 with --mpi 2.2.2.4.
# ---------------------------------------------------------------------------

# ---- EDIT: Grid Frontier tree + built executable ----
# Build: copy examples/Wilson_DW_spectrum.cc into the Grid tree, (re)run
# scripts/filelist so examples/Make.inc picks it up, then build in systems/Frontier;
# the binary lands in systems/Frontier/examples/.
GRID=$HOME/tmp/src/Grid_cleanedup_for_pullrequest        # Grid source root
EXE=${GRID}/systems/Frontier/examples/Wilson_DW_spectrum  # built executable

# ---- EDIT: which configs (glob, relative to the submit dir) ----
CONFIGS=( ckpoint_lat.* )       # e.g. ckpoint_lat.700 ckpoint_lat.710 ...

# ---- physics / solver parameters ----
MASS=0           # reference uses m = 0
LO=0.0; HI=2.0   # window Re(lambda) in [LO, HI]
NSHIFT=21        # real shifts tiling [LO,HI] (spacing 0.1)
KRYLOV=80        # RefinedArnoldi Krylov dim per shift (memory is ample on Frontier)
STOL=1e-11       # inner normal-equations CG tolerance (< ACCEPT)
SITER=30000      # inner CG max iterations
ACCEPT=1e-8      # keep modes with raw residual < ACCEPT
DEDUPE=1e-6      # merge modes within this distance in C

cd ${SLURM_SUBMIT_DIR}

cat << EOF > ./select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 4 5 6 7)
export NUMA_MAP=(3 3 1 1 0 0 2 2)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUM=\${NUMA_MAP[\$SLURM_LOCALID]}
export ROCR_VISIBLE_DEVICES=\$GPU
exec numactl -m \$NUM -N \$NUM \$*
EOF
chmod +x ./select_gpu

source ${GRID}/systems/Frontier/sourceme.sh
export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1
GRIDOPT="--accelerator-threads 8 --comms-overlap --shm 4096 --shm-mpi 0"

for CFG in "${CONFIGS[@]}"; do
    [ -e "${CFG}" ] || { echo "skip (missing): ${CFG}"; continue; }
    ID=$(echo "${CFG}" | grep -oE '[0-9]+' | tail -1)    # trailing number = config id / MD time
    OUT=evals_DW_${ID}.dat
    echo "=== ${CFG}  ->  ${OUT} ==="
    srun ./select_gpu ${EXE} \
        --grid 32.32.32.32 --mpi 2.2.2.4 ${GRIDOPT} \
        --config ${CFG} --mass ${MASS} \
        --window ${LO}:${HI} --nshift ${NSHIFT} --krylov ${KRYLOV} \
        --stol ${STOL} --siter ${SITER} --accept ${ACCEPT} --dedupe ${DEDUPE} \
        --tag ${ID} --out ${OUT}
done

echo "Done.  Per-config D_W eigenvalues (Re in [${LO},${HI}]) in evals_DW_*.dat"
