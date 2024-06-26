# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.29.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.29.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/v.reichi/projects/testing/cuda-seminar/2d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/v.reichi/projects/testing/cuda-seminar/2d/src

# Include any dependencies generated for this target.
include src/CMakeFiles/wavePacket.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/wavePacket.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/wavePacket.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/wavePacket.dir/flags.make

src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o: src/CMakeFiles/wavePacket.dir/flags.make
src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o: wavePacket.cpp
src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o: src/CMakeFiles/wavePacket.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/v.reichi/projects/testing/cuda-seminar/2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o"
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o -MF CMakeFiles/wavePacket.dir/wavePacket.cpp.o.d -o CMakeFiles/wavePacket.dir/wavePacket.cpp.o -c /Users/v.reichi/projects/testing/cuda-seminar/2d/src/wavePacket.cpp

src/CMakeFiles/wavePacket.dir/wavePacket.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/wavePacket.dir/wavePacket.cpp.i"
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/v.reichi/projects/testing/cuda-seminar/2d/src/wavePacket.cpp > CMakeFiles/wavePacket.dir/wavePacket.cpp.i

src/CMakeFiles/wavePacket.dir/wavePacket.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/wavePacket.dir/wavePacket.cpp.s"
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/v.reichi/projects/testing/cuda-seminar/2d/src/wavePacket.cpp -o CMakeFiles/wavePacket.dir/wavePacket.cpp.s

# Object files for target wavePacket
wavePacket_OBJECTS = \
"CMakeFiles/wavePacket.dir/wavePacket.cpp.o"

# External object files for target wavePacket
wavePacket_EXTERNAL_OBJECTS =

src/libwavePacket.a: src/CMakeFiles/wavePacket.dir/wavePacket.cpp.o
src/libwavePacket.a: src/CMakeFiles/wavePacket.dir/build.make
src/libwavePacket.a: src/CMakeFiles/wavePacket.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/v.reichi/projects/testing/cuda-seminar/2d/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libwavePacket.a"
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && $(CMAKE_COMMAND) -P CMakeFiles/wavePacket.dir/cmake_clean_target.cmake
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wavePacket.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/wavePacket.dir/build: src/libwavePacket.a
.PHONY : src/CMakeFiles/wavePacket.dir/build

src/CMakeFiles/wavePacket.dir/clean:
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src && $(CMAKE_COMMAND) -P CMakeFiles/wavePacket.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/wavePacket.dir/clean

src/CMakeFiles/wavePacket.dir/depend:
	cd /Users/v.reichi/projects/testing/cuda-seminar/2d/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/v.reichi/projects/testing/cuda-seminar/2d /Users/v.reichi/projects/testing/cuda-seminar/2d/src /Users/v.reichi/projects/testing/cuda-seminar/2d/src /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src /Users/v.reichi/projects/testing/cuda-seminar/2d/src/src/CMakeFiles/wavePacket.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/CMakeFiles/wavePacket.dir/depend

