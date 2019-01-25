#!/bin/bash
set -e -x

mkdir external
cd external

# Install slightly newer version of cmake 
curl -O https://cmake.org/files/v3.0/cmake-3.0.2.tar.gz
tar -zxvf cmake-3.0.2.tar.gz
cd cmake-3.0.2
./bootstrap
make -j2
make install
cd ..

# Install eigen3
curl -L -O http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
mv 3.2.10.tar.gz eigen-3.2.10.tar.gz
mkdir eigen-3.2.10
tar --directory=eigen-3.2.10 --strip-components=1 -zxvf eigen-3.2.10.tar.gz
cd eigen-3.2.10
mkdir build
cd build
cmake ..
make -j2 && make install
cd ..
cd ..

# Install GROMACS
curl -L -O http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz
tar -zxvf gromacs-2016.5.tar.gz
cd gromacs-2016.5

# Patch CMakeLists.txt
patch="if(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
patch+="    set(OLD_CMAKE_FIND_LIBRARY_SUFFIXES \${CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \".so\") \n"
patch+="    find_package(ZLIB QUIET) \n"
patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \${OLD_CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
patch+="else(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
patch+="    find_package(ZLIB QUIET) \n"
patch+="endif(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS)"
sed -i "s|find_package(ZLIB QUIET)|$patch|g" CMakeLists.txt

mkdir build
mkdir installed
export GMX_PATH=`pwd`/installed
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$GMX_PATH -DGMX_GPU=off -DBUILD_SHARED_LIBS=ON -DGMX_PREFER_STATIC_LIBS=ON -DCMAKE_CXX_FLAGS="-static-libstdc++" -DGMX_SIMD=SSE2
make -j2 && make install
cd ..
cd ..

# Out of external
cd ..


