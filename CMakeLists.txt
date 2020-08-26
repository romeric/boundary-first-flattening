cmake_minimum_required(VERSION 3.1.0)
project(bff)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# bff plus rectangle-bin-pack include paths
list(APPEND BFF_INCLUDE_DIRS "deps/rectangle-bin-pack/include")
list(APPEND BFF_INCLUDE_DIRS "mesh/include")
list(APPEND BFF_INCLUDE_DIRS "project/include")
include_directories(${BFF_INCLUDE_DIRS})

# Eigen include path
set(EIGEN_INCLUDE_DIR "$ENV{EIGEN_INCLUDE_DIR}" )
IF( NOT EIGEN_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Environment variable EIGEN_INCLUDE_DIR needs to be defined")
ENDIF()
include_directories("${EIGEN_INCLUDE_DIR}")

file(GLOB BFF_SOURCES "deps/rectangle-bin-pack/src/*.cpp" "linear-algebra/src/*.cpp" "mesh/src/*.cpp" "project/src/*.cpp")

if(MSVC)
    add_compile_options(
        /FS;
        /MD;
        )
endif()

# create bff static library
add_library(bff ${BFF_SOURCES})
target_link_libraries(bff ${MKL_LIBS})

if(MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "NDEBUG;"
    )
endif()

# build command line app
file(GLOB COMMAND_LINE_SOURCES "command-line/*.cpp")

add_executable(bff-command-line ${COMMAND_LINE_SOURCES})
target_link_libraries(bff-command-line bff)




