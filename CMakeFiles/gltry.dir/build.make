# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gero/gltry

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gero/gltry

# Include any dependencies generated for this target.
include CMakeFiles/gltry.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/gltry.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/gltry.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gltry.dir/flags.make

CMakeFiles/gltry.dir/src/main.cpp.o: CMakeFiles/gltry.dir/flags.make
CMakeFiles/gltry.dir/src/main.cpp.o: src/main.cpp
CMakeFiles/gltry.dir/src/main.cpp.o: CMakeFiles/gltry.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gero/gltry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/gltry.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gltry.dir/src/main.cpp.o -MF CMakeFiles/gltry.dir/src/main.cpp.o.d -o CMakeFiles/gltry.dir/src/main.cpp.o -c /home/gero/gltry/src/main.cpp

CMakeFiles/gltry.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltry.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gero/gltry/src/main.cpp > CMakeFiles/gltry.dir/src/main.cpp.i

CMakeFiles/gltry.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltry.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gero/gltry/src/main.cpp -o CMakeFiles/gltry.dir/src/main.cpp.s

CMakeFiles/gltry.dir/src/shaderprogram.cpp.o: CMakeFiles/gltry.dir/flags.make
CMakeFiles/gltry.dir/src/shaderprogram.cpp.o: src/shaderprogram.cpp
CMakeFiles/gltry.dir/src/shaderprogram.cpp.o: CMakeFiles/gltry.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gero/gltry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/gltry.dir/src/shaderprogram.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gltry.dir/src/shaderprogram.cpp.o -MF CMakeFiles/gltry.dir/src/shaderprogram.cpp.o.d -o CMakeFiles/gltry.dir/src/shaderprogram.cpp.o -c /home/gero/gltry/src/shaderprogram.cpp

CMakeFiles/gltry.dir/src/shaderprogram.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltry.dir/src/shaderprogram.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gero/gltry/src/shaderprogram.cpp > CMakeFiles/gltry.dir/src/shaderprogram.cpp.i

CMakeFiles/gltry.dir/src/shaderprogram.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltry.dir/src/shaderprogram.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gero/gltry/src/shaderprogram.cpp -o CMakeFiles/gltry.dir/src/shaderprogram.cpp.s

# Object files for target gltry
gltry_OBJECTS = \
"CMakeFiles/gltry.dir/src/main.cpp.o" \
"CMakeFiles/gltry.dir/src/shaderprogram.cpp.o"

# External object files for target gltry
gltry_EXTERNAL_OBJECTS =

gltry: CMakeFiles/gltry.dir/src/main.cpp.o
gltry: CMakeFiles/gltry.dir/src/shaderprogram.cpp.o
gltry: CMakeFiles/gltry.dir/build.make
gltry: /usr/lib/libglfw.so.3.3
gltry: /usr/lib/libGLEW.so
gltry: CMakeFiles/gltry.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gero/gltry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable gltry"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gltry.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gltry.dir/build: gltry
.PHONY : CMakeFiles/gltry.dir/build

CMakeFiles/gltry.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gltry.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gltry.dir/clean

CMakeFiles/gltry.dir/depend:
	cd /home/gero/gltry && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gero/gltry /home/gero/gltry /home/gero/gltry /home/gero/gltry /home/gero/gltry/CMakeFiles/gltry.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gltry.dir/depend

