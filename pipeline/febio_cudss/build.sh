#!/usr/bin/env bash
# Build the FEBioCuDSS plug-in (verified: docs/FINDINGS.md sections 28, 30). Needs the CUDA toolkit
# (nvcc/CUDAToolkit), a FEBio build (pipeline/build_febio.sh) and cuDSS.
#
# usage: FEBIO=$HOME/febio/FEBio CUDSS=$HOME/.pylibs_cudss/nvidia/cu12 bash pipeline/febio_cudss/build.sh
# If the toolkit is newer than CUDA 12 (cuDSS wheel is cu12), also set
#   CUDART_LIB=.../nvidia/cuda_runtime/lib/libcudart.so.12  CUDART_INCLUDE=.../nvidia/cuda_runtime/include
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
FEBIO=${FEBIO:-$HOME/febio/FEBio}
CUDSS=${CUDSS:-$HOME/.pylibs_cudss/nvidia/cu12}
# the pip wheel ships only libcudss.so.0; the linker wants the unversioned name
[ -e "$CUDSS/lib/libcudss.so" ] || ln -s libcudss.so.0 "$CUDSS/lib/libcudss.so"
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release \
      -DFEBIO_ROOT="$FEBIO" -DFEBIO_LIB_DIR="$FEBIO/build/lib" -DCUDSS_ROOT="$CUDSS" ${CUDART_LIB:+-DCUDART_LIB="$CUDART_LIB"} ${CUDART_INCLUDE:+-DCUDART_INCLUDE="$CUDART_INCLUDE"}
cmake --build "$HERE/build" -j "$(nproc)"
echo "plug-in: $HERE/build/libFEBioCuDSS.so"
echo "run:     LD_LIBRARY_PATH=$CUDSS/lib:\$LD_LIBRARY_PATH febio4 -import $HERE/build/libFEBioCuDSS.so -i model.feb"
