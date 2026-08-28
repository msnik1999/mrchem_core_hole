#!/bin/bash -l

if [ -z "$1" ]; then
    # If not, prompt the user for input
    echo "This script requires an NRIS account for compilation on a compute node"
    read -p "Please state your account (nnxxxxk): " account
else
    # Otherwise, use the argument
    account="$1"
fi

mrchem_dir="$(pwd)"
source ${mrchem_dir}/tools/olivia.env

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

./setup --prefix=${install_dir} --omp --mpi ${build_dir} && \
cd ${build_dir} && \
srun --cpus-per-task=4 --mem-per-cpu=2G --time=01:00:00 --account=$account bash -c \
"source ${mrchem_dir}/tools/olivia.env ; \
make -j4 ; \
OMP_NUM_THREADS=4 ctest -L unit --output-on-failure ; \
make install -j4"

exit 0
