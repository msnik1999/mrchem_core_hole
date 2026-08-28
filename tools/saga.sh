#!/bin/bash -l

if [ `hostname | grep -i login | wc -l` == 0 ]; then
   echo "This script MUST be run on the LOGIN node as it requires internet access"
   exit 1
fi

mrchem_dir="$(pwd)"
source ${mrchem_dir}/tools/saga.env

cd ${mrchem_dir}
build_dir=${mrchem_dir}/build
install_dir=${mrchem_dir}/install

if [ -d "${build_dir}" ]; then
    echo "Build directory already exists."
    read -p $'\nWould you like to delete the existing build and install folders and continue? (y/N): ' response
    case "${response,}" in
        y|yes)
            echo "Removing old build and install directories..."
            rm -rf "${build_dir}"
            rm -rf "${install_dir}"
            ;;
        *)
            echo "Aborting build to prevent overwriting..."
            exit 1
            ;;
    esac
fi

./setup --prefix=${install_dir} --omp --mpi --cxx=mpicxx ${build_dir} && \
cd ${build_dir} && \
make && \
OMP_NUM_THREADS=1 ctest -L unit --output-on-failure && \
make install

exit 0
