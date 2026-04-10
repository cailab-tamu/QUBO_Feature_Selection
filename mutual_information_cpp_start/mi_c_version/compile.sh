#!/bin/bash

# Compiler
CXX=g++
MPICXX=mpic++

# Flags
CXXFLAGS_RELEASE="-O2 -Wall -I/usr/include/eigen3 -lm"
CXXFLAGS_DEBUG="-g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I/usr/include/eigen3 -lm"
OPENMPFLAGS="-fopenmp"
MPIFLAGS=""

# Build type (default is Release)
if [ "$DEBUG" == "yes" ]; then
    CXXFLAGS=$CXXFLAGS_DEBUG
else
    CXXFLAGS=$CXXFLAGS_RELEASE
fi

# Source Files
SERIAL_SRC="mi_serial.cpp"
OPENMP_SRC="mi_openmp.cpp"
MPI_SRC="mi_mpi.cpp"
MMIO_SRC="mmio.c"
MMIO_HEADER="mmio.h"

# Targets
SERIAL_TARGET="mi_serial"
OPENMP_TARGET="mi_openmp"
MPI_TARGET="mi_mpi"

# Compile serial version
$CXX $CXXFLAGS -o $SERIAL_TARGET $SERIAL_SRC $MMIO_SRC
echo "Compiled $SERIAL_TARGET"

# Compile OpenMP version
$CXX $CXXFLAGS $OPENMPFLAGS -o $OPENMP_TARGET $OPENMP_SRC $MMIO_SRC
echo "Compiled $OPENMP_TARGET"

# Compile MPI version
$MPICXX $CXXFLAGS $MPIFLAGS -o $MPI_TARGET $MPI_SRC $MMIO_SRC
echo "Compiled $MPI_TARGET"

# Clean targets
clean() {
    rm -f $SERIAL_TARGET $OPENMP_TARGET $MPI_TARGET
    echo "Cleaned up compiled targets"
}

# If argument is "clean", remove compiled targets
if [ "$1" == "clean" ]; then
    clean
fi

