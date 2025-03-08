#!/bin/bash
set -e -x

# assume that conda is installed and venv environment is created
conda activate ./venv

CWD=`pwd` # Current working directory

# Build and Install GROMACS
cd external
if [ -d gmx_installed ]; then
    rm -rf gmx_installed
fi
mkdir gmx_installed
cd gromacs

if [ -d build ]; then
    rm -rf build
fi

# add conda library to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$CWD/venv/lib

mkdir build && cd build
cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack \
        -DCMAKE_INSTALL_PREFIX=${CWD}/external/gmx_installed
make -j 10
make install

# Build and install gmx_clusterByFeatures in in-place mode
cd ${CWD}
export GMX_SRC=${CWD}/external/gromacs
export GMX_INSTALL=${CWD}/external/gmx_installed
GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip install -ve .
