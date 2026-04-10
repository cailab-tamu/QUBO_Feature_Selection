CXX = g++

# Default SYS
SYS ?= linux


LOCAL_EIGEN_WIN = C:/Users/ssromerogon/Documents/vscode_working_dir/eigen # Local installation
LOCAL_EIGEN_LIN = /usr/include/eigen3

# Flags for Windows
CXXFLAGS_WINDOWS_RELEASE = -O2 -Wall -I $(LOCAL_EIGEN_WIN) -lm
CXXFLAGS_WINDOWS_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I $(LOCAL_EIGEN_WIN) -lm

# Flags for Linux
CXXFLAGS_LINUX_RELEASE = -O2 -Wall -I $(LOCAL_EIGEN_LIN) -lm
CXXFLAGS_LINUX_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I $(LOCAL_EIGEN_LIN) -lm

# OpenMP flags
OPENMPFLAGS = -fopenmp

# Select SYS-specific flags and print platform info
ifeq ($(SYS),windows)
    ifeq ($(DEBUG),yes)
        CXXFLAGS = $(CXXFLAGS_WINDOWS_DEBUG)
    else
        CXXFLAGS = $(CXXFLAGS_WINDOWS_RELEASE)
    endif
    PLATFORM_MESSAGE = "Building on Windows"
    EIGEN_PATH = $(LOCAL_EIGEN_WIN)
else ifeq ($(SYS),linux)
    ifeq ($(DEBUG),yes)
        CXXFLAGS = $(CXXFLAGS_LINUX_DEBUG)
    else
        CXXFLAGS = $(CXXFLAGS_LINUX_RELEASE)
    endif
    PLATFORM_MESSAGE = "Building on Linux"
    EIGEN_PATH = $(LOCAL_EIGEN_LIN)
else
    $(error Unsupported platform: $(SYS))
endif

# Targets
TARGETS = mi_serial mi_openmp

# Rules
all: $(TARGETS)
	@echo "$(PLATFORM_MESSAGE)"
	@echo "Eigen Path: $(EIGEN_PATH)"
	@echo "CXXFLAGS: $(CXXFLAGS)"


mi_serial: mi_serial.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) -o $@ mi_serial.cpp mmio.c

mi_openmp: mi_openmp.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ mi_openmp.cpp mmio.c

clean:
	rm -f $(TARGETS)

.PHONY: all clean  # Add .PHONY to prevent issues with files named 'all' or 'clean'
