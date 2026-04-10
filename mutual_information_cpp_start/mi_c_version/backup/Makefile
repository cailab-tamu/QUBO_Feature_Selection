# Compiler
CXX = g++

# Default SYS
SYS ?= linux
LOCAL_EIGEN_WIN = C:/Users/ssromerogon/Documents/vscode_working_dir/eigen # Pulled repo
#LOCAL_EIGEN_WIN = C:/Users/ssromerogon/Documents/vscode_working_dir/scGEAToolbox_addson/QUBO_FS_python/mi_c_version/eigen3/include/eigen3 # Local installation
LOCAL_EIGEN_LIN = /usr/include/eigen3

# Flags for Windows
CXXFLAGS_WINDOWS_RELEASE = -O2 -Wall -I $(LOCAL_EIGEN_WIN) -lm
CXXFLAGS_WINDOWS_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I $(LOCAL_EIGEN_WIN) -lm

# Flags for Linux
CXXFLAGS_LINUX_RELEASE = -O2 -Wall -I $(LOCAL_EIGEN_LIN) -lm
CXXFLAGS_LINUX_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I $(LOCAL_EIGEN_LIN) -lm

# OpenMP flags
OPENMPFLAGS = -fopenmp

# Select SYS-specific flags
ifeq ($(SYS),windows)
    ifeq ($(DEBUG),yes)
        CXXFLAGS = $(CXXFLAGS_WINDOWS_DEBUG)
    else
        CXXFLAGS = $(CXXFLAGS_WINDOWS_RELEASE)
    endif
else ifeq ($(SYS),linux)
    ifeq ($(DEBUG),yes)
        CXXFLAGS = $(CXXFLAGS_LINUX_DEBUG)
    else
        CXXFLAGS = $(CXXFLAGS_LINUX_RELEASE)
    endif
else
    $(error Unsupported platform: $(SYS))
endif

# Targets
TARGETS = mi_serial mi_openmp

# Rules
all: $(TARGETS)

mi_serial: mi_serial.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) -o $@ mi_serial.cpp mmio.c

mi_openmp: mi_openmp.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ mi_openmp.cpp mmio.c

clean:
	rm -f $(TARGETS)
