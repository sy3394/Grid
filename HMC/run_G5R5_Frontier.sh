#!/bin/bash -l
#SBATCH --job-name=DWF_G5R5
#SBATCH --partition=extended
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=10:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --out=out/G5R5_DWF-%j.out

CPU_BIND="map_ldom:3,3,1,1,0,0,2,2"
MEM_BIND="map_mem:3,3,1,1,0,0,2,2"
echo $CPU_BIND

cat << 'EOF' > ./select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 4 5 6 7)
export NUMA_MAP=(3 3 1 1 0 0 2 2)
export GPU=${GPU_MAP[$SLURM_LOCALID]}
export NUM=${NUMA_MAP[$SLURM_LOCALID]}
export ROCR_VISIBLE_DEVICES=$GPU
echo RANK $SLURM_LOCALID using GPU $GPU
echo NUMA $SLURM_LOCALID using NUMA ${NUM}
exec numactl -m $NUM -N $NUM $*
EOF
chmod +x ./select_gpu

root=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier
BINARY=${root}/HMC/Compute_DWF_G5R5
source ${root}/sourceme.sh
module list

export OMP_NUM_THREADS=7
export MPICH_SMP_SINGLE_COPY_MODE=CMA
export MPICH_GPU_SUPPORT_ENABLED=1

### Prepare input / output directories
traj=702
n_tot=1
fpath=$(pwd)
fname=ckpoint_EODWF_lat_smr
outpath=$(pwd)/eigen

for (( i=traj; i<traj+n_tot; i++ )); do
    mkdir -p $outpath/$i
done

# Compute_DWF_G5R5 reads LanParams.xml from the current working directory (hardcoded).
cat > LanParams.xml << EOF
<?xml version="1.0"?>
<grid>
  <LanczosParameters>
    <mass>1.0</mass>
    <M5>1.8</M5>
    <Ls>48</Ls>
    <Nstop>16</Nstop>
    <Nk>16</Nk>
    <Np>32</Np>
    <ChebyLow>0.003</ChebyLow>
    <ChebyHigh>68</ChebyHigh>
    <ChebyOrder>201</ChebyOrder>
    <StartTrajectory>${traj}</StartTrajectory>
    <Trajectories>${n_tot}</Trajectories>
    <fpath>${fpath}</fpath>
    <fname>${fname}</fname>
    <outpath>${outpath}</outpath>
  </LanczosParameters>
  <WilsonFlow>
    <is_flow>false</is_flow>
    <take_meas>false</take_meas>
    <steps>0</steps>
    <step_size>0.01</step_size>
    <meas_interval>400</meas_interval>
    <maxTau>0</maxTau>
    <path>${fpath}/dnsty</path>
  </WilsonFlow>
</grid>
EOF
cat LanParams.xml

echo ======================
echo Running G5R5 eigenvalue + density
echo ======================

vol=32.32.32.32
mpi=4.2.2.2
PARAMS="--mpi $mpi --accelerator-threads 32 --shm 2048 --shm-mpi 0 --grid $vol"
srun ./select_gpu $BINARY $PARAMS
