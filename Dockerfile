FROM quay.io/pypa/manylinux_2_28_x86_64

COPY gmx_clusterByFeatures /app-src/gmx_clusterByFeatures
COPY external/gromacs /app-src/external/gromacs
COPY src /app-src/src
COPY .git /app-src/.git
COPY pyCode2Hex.py /app-src/pyCode2Hex.py
COPY setup.py /app-src/setup.py

WORKDIR /app-src/external
RUN mkdir gmx_installed 
WORKDIR /app-src/external/gromacs
RUN rm -rf build 
RUN mkdir build
WORKDIR /app-src/external/gromacs/build

ENV GMX_PATH=/app-src/external/gmx_installed
ENV GMX_SRC=/app-src/external/gromacs

RUN cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_PATH}
RUN make
RUN make install

COPY scripts /app-src/scripts
WORKDIR /app-src
RUN bash -i scripts/build_wheels.sh

CMD ["echo", "testing"]