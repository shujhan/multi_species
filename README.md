# multi_species
module load nvhpc-openmpi4/22.3 cuda/11.3

module load eigen boost

cmake .. -D CMAKE_CXX_COMPILER=nvc++

make

After getting the compiled file farrsight_cpu and farrsight_gpu, use python interface and .sh file to run different jobs.
