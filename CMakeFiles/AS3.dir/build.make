# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3

# Include any dependencies generated for this target.
include CMakeFiles/AS3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/AS3.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/AS3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/AS3.dir/flags.make

CMakeFiles/AS3.dir/main.cpp.o: CMakeFiles/AS3.dir/flags.make
CMakeFiles/AS3.dir/main.cpp.o: main.cpp
CMakeFiles/AS3.dir/main.cpp.o: CMakeFiles/AS3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/AS3.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/AS3.dir/main.cpp.o -MF CMakeFiles/AS3.dir/main.cpp.o.d -o CMakeFiles/AS3.dir/main.cpp.o -c /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/main.cpp

CMakeFiles/AS3.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AS3.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/main.cpp > CMakeFiles/AS3.dir/main.cpp.i

CMakeFiles/AS3.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AS3.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/main.cpp -o CMakeFiles/AS3.dir/main.cpp.s

# Object files for target AS3
AS3_OBJECTS = \
"CMakeFiles/AS3.dir/main.cpp.o"

# External object files for target AS3
AS3_EXTERNAL_OBJECTS =

AS3: CMakeFiles/AS3.dir/main.cpp.o
AS3: CMakeFiles/AS3.dir/build.make
AS3: /usr/local/lib/libglfw3.a
AS3: CMakeFiles/AS3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable AS3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AS3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/AS3.dir/build: AS3
.PHONY : CMakeFiles/AS3.dir/build

CMakeFiles/AS3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/AS3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/AS3.dir/clean

CMakeFiles/AS3.dir/depend:
	cd /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3 /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3 /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3 /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3 /Users/filippocasari/Dropbox/Mac/Documents/ParticleMethods/Assignment3/CMakeFiles/AS3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/AS3.dir/depend
