# FEBioCuDSS plug-in. See README.md and docs/FINDINGS.md section 30.
cmake_minimum_required(VERSION 3.18)
project(FEBioCuDSS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(CUDAToolkit REQUIRED)
set(FEBIO_ROOT "" CACHE PATH "FEBio source root (the directory containing FECore/)")
set(FEBIO_LIB_DIR "" CACHE PATH "FEBio build library directory (contains libfecore)")
set(CUDSS_ROOT "" CACHE PATH "cuDSS root (contains include/cudss.h and lib/libcudss.so)")

add_library(FEBioCuDSS SHARED FEBioCuDSS.cpp)
target_compile_definitions(FEBioCuDSS PRIVATE LINUX)
set(CUDART_INCLUDE "" CACHE PATH "CUDA runtime headers matching CUDART_LIB (e.g. nvidia/cuda_runtime/include)")
target_include_directories(FEBioCuDSS PRIVATE ${FEBIO_ROOT} ${CUDSS_ROOT}/include ${CUDART_INCLUDE} ${CUDAToolkit_INCLUDE_DIRS})
target_link_directories(FEBioCuDSS PRIVATE ${FEBIO_LIB_DIR} ${CUDSS_ROOT}/lib)
# cuDSS from the cu12 wheel needs the CUDA 12 runtime; on a host whose toolkit is newer,
# point CUDART_LIB at a libcudart.so.12 (e.g. the one torch ships) so one runtime is loaded.
set(CUDART_LIB "" CACHE FILEPATH "explicit libcudart to link (default: the toolkit's)")
if(CUDART_LIB)
  target_link_libraries(FEBioCuDSS PRIVATE fecore cudss ${CUDART_LIB})
else()
  target_link_libraries(FEBioCuDSS PRIVATE fecore cudss CUDA::cudart)
endif()
