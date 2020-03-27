# This file specifies how the project should be built, using CMake

##################################################################################
# Configurations
# Set project name, version and language
cmake_minimum_required(VERSION 3.15)
project(cycleCount VERSION 1.00 LANGUAGES CXX)

##################################################################################
# Files location
# All .cpp files in src/
set(SOURCES
        src/cycleCount.cpp
        src/eigenvalues.cpp
        )
# Default name for the library built from src/*.cpp
set(LIBRARY_NAME cycleCount)

##################################################################################
# Build
# Compile all sources into a library.
add_library(${LIBRARY_NAME} OBJECT ${SOURCES})
# Lib needs its header files, and users of the library must also see these.
target_include_directories(${LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/header)
# Set compiler options you want
set(CMAKE_CXX_FLAGS "--coverage")

##################################################################################
# Executable
# Name of the executable and location of file
add_executable(main app/main.cpp)
# Link the executable to library
target_link_libraries(main PRIVATE ${LIBRARY_NAME})

set_target_properties(
        ${LIBRARY_NAME} main
        PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

