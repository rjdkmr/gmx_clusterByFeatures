#!/bin/bash
set -e -x

gcc --version

yum -y install epel-release
yum -y install cairo
yum -y install cairo-devel
yum -y install libxml2-devel
yum -y remove cmake28
yum -y install fftw3 fftw3-devel
yum -y install openblas openblas-devel blas blas-devel atlas atlas-devel lapack lapack-devel

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

# Compile wheels
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${GMX_PATH}/lib64/pkgconfig/"
for PYBIN in /opt/python/cp3*/bin; do
    "${PYBIN}/pip" install -r /io/dev-requirements.txt
    "${PYBIN}/pip" wheel -w /io/wheelhouse/ --no-deps --no-cache-dir /io/
done

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse/*.whl; do
    auditwheel show "$whl" 
    auditwheel repair "$whl" -w /io/wheelhouse/
done


# Install packages and test
for PYBIN in /opt/python/cp3*/bin/; do
    "${PYBIN}pip" install gmx_clusterByFeatures --no-index -f /io/wheelhouse
    "${PYBIN}python" -c "import gmx_clusterByFeatures; print('=====\nTEST -- gmx_clusterByFeatures GROMACS version: ', gmx_clusterByFeatures.gmx_version, '\n=====')"
done

