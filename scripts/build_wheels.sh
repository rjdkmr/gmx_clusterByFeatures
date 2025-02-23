export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${GMX_INSATLL}/lib64/pkgconfig/"

cd /

for PYBIN in /opt/python/cp3*/bin/; do
    echo "PYBIN: ${PYBIN}"
done

export GMX_INSTALL=/app-src/external/gmx_installed
export GMX_SRC=/app-src/external/gromacs

pyTags=("cp38-cp38" "cp39-cp39" "cp310-cp310" "cp311-cp311")
for pyTag in ${pyTags[@]}; do
    PYBIN="/opt/python/${pyTag}/bin"
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} "${PYBIN}/python" -m pip install -r /io/dev-requirements.txt
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} "${PYBIN}/python" -m pip wheel -w /io/wheelhouse/ --no-deps --no-cache-dir /io
done

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse/*.whl; do
    auditwheel show "$whl"
    auditwheel repair "$whl" -w /io/wheelhouse/
done


# Install packages and test
for pyTag in ${pyTags[@]}; do
    PYBIN="/opt/python/${pyTag}/bin"
    "${PYBIN}/python" -m pip install gmx_clusterByFeatures --no-index -f /io/wheelhouse
    "${PYBIN}/python" -c "import gmx_clusterByFeatures; print('=====\nTEST -- gmx_clusterByFeatures GROMACS version: ', gmx_clusterByFeatures.gmx_version, '\n=====')"
done
