cmake_minimum_required(VERSION 3.1)
project(multibody)
add_definitions(-DHAVE_CBLAS=1)

if (WIN32)
    find_package(OpenBLAS REQUIRED)
    set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()


find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)

add_executable(multibody Dyn_simulator.cpp)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS)

message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})

target_link_libraries(multibody ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} xtensor xtensor::optimize xtensor::use_xsimd xtensor-blas)

#${PYTHON_LIBRARIES}         Python3::NumPy)

