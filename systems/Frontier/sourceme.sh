. /autofs/nccs-svm1_home1/paboyle/Crusher/Grid/spack/share/spack/setup-env.sh
spack load c-lime gmp@6.2.1
module load emacs 
#module load PrgEnv-gnu
module load PrgEnv-cray-amd/8.5.0
module load rocm #/6.2.0
module load cray-mpich
#module load gmp-6.2.1-gcc-7.5.0-h5wpx5n  #gmp
module load cray-fftw
module load craype-accel-amd-gfx90a
export LD_LIBRARY_PATH=/opt/gcc/mpfr/3.1.4/lib:$LD_LIBRARY_PATH
#Hack for lib
#export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
