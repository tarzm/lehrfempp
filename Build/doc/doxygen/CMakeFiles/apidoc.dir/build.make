# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/magina/Documents/lehrfempp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/magina/Documents/lehrfempp/Build

# Utility rule file for apidoc.

# Include the progress variables for this target.
include doc/doxygen/CMakeFiles/apidoc.dir/progress.make

doc/doxygen/CMakeFiles/apidoc: documentation


documentation: ../doc/doxygen/Doxyfile
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/u/magina/Documents/lehrfempp/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating ../../documentation"
	cd /u/magina/Documents/lehrfempp/doc/doxygen && /usr/bin/cmake -E echo_append Building\ API\ Documentation...
	cd /u/magina/Documents/lehrfempp/doc/doxygen && /usr/bin/doxygen Doxyfile
	cd /u/magina/Documents/lehrfempp/doc/doxygen && /usr/bin/cmake -E echo Done.

apidoc: doc/doxygen/CMakeFiles/apidoc
apidoc: documentation
apidoc: doc/doxygen/CMakeFiles/apidoc.dir/build.make

.PHONY : apidoc

# Rule to build all files generated by this target.
doc/doxygen/CMakeFiles/apidoc.dir/build: apidoc

.PHONY : doc/doxygen/CMakeFiles/apidoc.dir/build

doc/doxygen/CMakeFiles/apidoc.dir/clean:
	cd /u/magina/Documents/lehrfempp/Build/doc/doxygen && $(CMAKE_COMMAND) -P CMakeFiles/apidoc.dir/cmake_clean.cmake
.PHONY : doc/doxygen/CMakeFiles/apidoc.dir/clean

doc/doxygen/CMakeFiles/apidoc.dir/depend:
	cd /u/magina/Documents/lehrfempp/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/magina/Documents/lehrfempp /u/magina/Documents/lehrfempp/doc/doxygen /u/magina/Documents/lehrfempp/Build /u/magina/Documents/lehrfempp/Build/doc/doxygen /u/magina/Documents/lehrfempp/Build/doc/doxygen/CMakeFiles/apidoc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/doxygen/CMakeFiles/apidoc.dir/depend

