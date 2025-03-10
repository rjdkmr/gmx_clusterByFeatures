#!/bin/bash
set -e -x

CWD=`pwd`
cd external
mkdir gmx_installed 
cd gromacs

ls -lrta

if [ -d build ]; then
    rm -rf build
fi

mkdir build
cd build

export GMX_INSTALL=${CWD}/external/gmx_installed
export GMX_SRC=${CWD}/external/gromacs

cmake -DGMX_SIMD=NONE -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_INSTALL} ..
make
make install

cd $CWD
