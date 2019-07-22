cmake_minimum_required(VERSION 3.9)
project(specCoarsen)

# EIGEN
# ---------------------------------------------------------------------------------
add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE external/eigen)
if(MSVC)
	target_sources(eigen INTERFACE ${CMAKE_CURRENT_LIST_DIR}/external/eigen/debug/msvc/eigen.natvis)
endif()

# libigl
# ---------------------------------------------------------------------------------
add_library(libigl INTERFACE)
target_include_directories(libigl INTERFACE external/libigl/include)

# Spectra
# ---------------------------------------------------------------------------------
add_library(spectra INTERFACE)
target_include_directories(spectra INTERFACE external/spectra/include)

# Spectral coarsening
# =================================================================================
file(GLOB SRCFILES main.cpp src/*.cpp src/*.h)
add_executable(specCoarsen ${SRCFILES})
target_link_libraries(specCoarsen
	PUBLIC eigen
	PUBLIC libigl
	PUBLIC spectra)
