#!/bin/bash
set -e -x

CWD=`pwd`

brew install gcc@14
brew install libomp
brew install gsl
brew install fftw
brew install pyenv
brew install eigen
brew cleanup

cd external
mkdir gmx_installed 
cd gromacs

if [ -d build ]; then
    rm -rf build
fi

mkdir build
cd build

export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"

export GMX_INSTALL=${CWD}/external/gmx_installed
export GMX_SRC=${CWD}/external/gromacs

export CC=gcc-14
export CXX=g++-14
cmake -DCMAKE_CC_COMIPLER=gcc-14 -DCMAKE_CXX_COMIPLER=g++-14 -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_INSTALL} ..
make
make install

cd $CWD

eval "$(pyenv init -)"
pyenv install --list
PYVERS=("3.9" "3.10" "3.11" "3.12")
PYTHONS=()
for PYVER in ${PYVERS[@]}
do
    # Install the latest release of the specified Python version using pyenv.
    PYVER="$(pyenv install --list | grep -E "^\\s*$PYVER" | sort -n -t. -k3 | tail -n1)"
    pyenv install $PYVER
    pyenv global $PYVER
    PYTHONS+=($PYVER)
    echo $PYVER
    echo $(python --version)
    python -m pip install --upgrade pip
    python -m python -m pip install -r {project}/dev-requirements.txt
    python -m pip install delocate
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip install -v --no-deps --no-cache-dir .
    install_name_tool -change @rpath/libgromacs.2.dylib ${GMX_INSTALL}/lib/libgromacs.dylib build/lib.*/gmx_clusterByFeatures/*.so
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip wheel -v -w wheels/ --no-deps --no-cache-dir .
done

delocate-listdeps wheels/*.whl
delocate-wheel -w fixed_wheels -v wheels/*.whl
delocate-listdeps fixed_wheels/*.whl

ls -lrt wheels/

cd $CWD