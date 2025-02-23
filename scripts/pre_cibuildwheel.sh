#!/bin/bash
set -e -x

CWD=`pwd`
cd external
mkdir gmx_installed 
cd gromacs

if [ -d build ]; then
    rm -rf build
fi

mkdir build
cd buiild

export GMX_INSTALL=/app-src/external/gmx_installed
export GMX_SRC=/app-src/external/gromacs

cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_INSTALL}
make
make install

cd $CWD
