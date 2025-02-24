#!/bin/bash
set -e -x

brew cask uninstall oclint
brew install gcc@9
brew install gsl
brew install fftw
brew outdated pyenv || brew upgrade pyenv
brew upgrade cmake
brew install eigen
brew cleanup

CWD=`pwd`
cd external
mkdir gmx_installed 
cd gromacs

if [ -d build ]; then
    rm -rf build
fi

mkdir build
cd build

export GMX_INSTALL=${CWD}/external/gmx_installed
export GMX_SRC=${CWD}/external/gromacs

export C=gcc
export CXX=g++
cmake -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_INSTALL} ..
make
make install

cd /project

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
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip build -v --no-deps --no-cache-dir .
    install_name_tool -change @rpath/libgromacs.2.dylib ${GMX_INSTALL}/lib/libgromacs.dylib build/lib.*/gmx_clusterByFeatures/*.so
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip wheel -v -w wheels/ --no-deps --no-cache-dir .
done

delocate-listdeps wheels/*.whl
delocate-wheel -w fixed_wheels -v wheels/*.whl
delocate-listdeps fixed_wheels/*.whl

ls -lrt wheels/

cd $CWD