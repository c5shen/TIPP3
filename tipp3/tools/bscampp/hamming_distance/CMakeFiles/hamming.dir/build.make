# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /home/ellie/anaconda3/lib/python3.8/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/ellie/anaconda3/lib/python3.8/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/ellie/easystore/BATCH-SCAMPP/scripts

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/ellie/easystore/BATCH-SCAMPP/scripts

# Include any dependencies generated for this target.
include CMakeFiles/hamming.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/hamming.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/hamming.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/hamming.dir/flags.make

CMakeFiles/hamming.dir/new_hamming.cpp.o: CMakeFiles/hamming.dir/flags.make
CMakeFiles/hamming.dir/new_hamming.cpp.o: new_hamming.cpp
CMakeFiles/hamming.dir/new_hamming.cpp.o: CMakeFiles/hamming.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ellie/easystore/BATCH-SCAMPP/scripts/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/hamming.dir/new_hamming.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/hamming.dir/new_hamming.cpp.o -MF CMakeFiles/hamming.dir/new_hamming.cpp.o.d -o CMakeFiles/hamming.dir/new_hamming.cpp.o -c /media/ellie/easystore/BATCH-SCAMPP/scripts/new_hamming.cpp

CMakeFiles/hamming.dir/new_hamming.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hamming.dir/new_hamming.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/ellie/easystore/BATCH-SCAMPP/scripts/new_hamming.cpp > CMakeFiles/hamming.dir/new_hamming.cpp.i

CMakeFiles/hamming.dir/new_hamming.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hamming.dir/new_hamming.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/ellie/easystore/BATCH-SCAMPP/scripts/new_hamming.cpp -o CMakeFiles/hamming.dir/new_hamming.cpp.s

# Object files for target hamming
hamming_OBJECTS = \
"CMakeFiles/hamming.dir/new_hamming.cpp.o"

# External object files for target hamming
hamming_EXTERNAL_OBJECTS =

hamming: CMakeFiles/hamming.dir/new_hamming.cpp.o
hamming: CMakeFiles/hamming.dir/build.make
hamming: /usr/lib/gcc/x86_64-linux-gnu/10/libgomp.so
hamming: /usr/lib/x86_64-linux-gnu/libpthread.so
hamming: CMakeFiles/hamming.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/ellie/easystore/BATCH-SCAMPP/scripts/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable hamming"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hamming.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/hamming.dir/build: hamming
.PHONY : CMakeFiles/hamming.dir/build

CMakeFiles/hamming.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/hamming.dir/cmake_clean.cmake
.PHONY : CMakeFiles/hamming.dir/clean

CMakeFiles/hamming.dir/depend:
	cd /media/ellie/easystore/BATCH-SCAMPP/scripts && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/ellie/easystore/BATCH-SCAMPP/scripts /media/ellie/easystore/BATCH-SCAMPP/scripts /media/ellie/easystore/BATCH-SCAMPP/scripts /media/ellie/easystore/BATCH-SCAMPP/scripts /media/ellie/easystore/BATCH-SCAMPP/scripts/CMakeFiles/hamming.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/hamming.dir/depend
