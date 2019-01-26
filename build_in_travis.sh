#!/bin/bash
set -e -x

mkdir external
cd external

# Install GROMACS
curl -L -O http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz
tar -zxf gromacs-2016.5.tar.gz
cd gromacs-2016.5

# Patch CMakeLists.txt
if [[ $1 != 'osx' ]]; then
  patch="if(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
  patch+="    set(OLD_CMAKE_FIND_LIBRARY_SUFFIXES \${CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
  patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \".so\") \n"
  patch+="    find_package(ZLIB QUIET) \n"
  patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \${OLD_CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
  patch+="else(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
  patch+="    find_package(ZLIB QUIET) \n"
  patch+="endif(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS)"
  sed -i "s|find_package(ZLIB QUIET)|$patch|g" CMakeLists.txt
fi

mkdir build
cd build
if [[ $1 != 'osx' ]]; then
  cmake .. -DGMX_GPU=off -DBUILD_SHARED_LIBS=ON -DGMX_PREFER_STATIC_LIBS=ON -DCMAKE_CXX_FLAGS="-static-libstdc++" -DGMX_SIMD=SSE2
else
  cmake .. -DGMX_GPU=off -DBUILD_SHARED_LIBS=ON -DGMX_PREFER_STATIC_LIBS=ON -DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -lomp" -DGMX_SIMD=SSE2
fi
make -j2 && sudo make install
cd ..
cd ..

# Out of external
cd ..


