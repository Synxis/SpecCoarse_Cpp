cmake_minimum_required(VERSION 3.9)
project(specCoarsen)

# OPENGL
# ---------------------------------------------------------------------------------
find_package(OpenGL)

# GLAD
# ---------------------------------------------------------------------------------
add_subdirectory(external/glad)

# GLFW
# ---------------------------------------------------------------------------------
option(GLFW_BUILD_EXAMPLES "" OFF)
option(GLFW_BUILD_TESTS    "" OFF)
option(GLFW_BUILD_DOCS     "" OFF)
option(GLFW_INSTALL        "" OFF)
option(GLFW_VULKAN_STATIC  "" OFF)
option(GLFW_USE_OSMESA     "" OFF)
option(GLFW_USE_HYBRID_HPG "" OFF)
option(GLFW_USE_WAYLAND    "" OFF)
option(GLFW_USE_MIR        "" OFF)
add_subdirectory(external/glfw)

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
target_include_directories(specCoarsen PUBLIC src)
target_link_libraries(specCoarsen
	PUBLIC eigen
	PUBLIC libigl
	PUBLIC spectra
	PUBLIC OpenGL::GL
	PUBLIC glad
	PUBLIC glfw)
