# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /cluster/apps/gcc-9.3.0/cmake-3.20.3-qrkwufs7msvgysj7lm64kcmd2hhj37m4/bin/cmake

# The command to remove a file.
RM = /cluster/apps/gcc-9.3.0/cmake-3.20.3-qrkwufs7msvgysj7lm64kcmd2hhj37m4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cluster/home/tamaurer/lehrfempp/lehrfempp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cluster/home/tamaurer/lehrfempp/lehrfempp/build

# Include any dependencies generated for this target.
include projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/compiler_depend.make

# Include the progress variables for this target.
include projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/progress.make

# Include the compile flags for this target's objects.
include projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/flags.make

projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o: projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/flags.make
projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o: ../projects/dgfem/fe_experiments.cc
projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o: projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cluster/home/tamaurer/lehrfempp/lehrfempp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o"
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem && /cluster/apps/gcc-4.8.5/gcc-9.3.0-mjmsvu3362usxjr7tas4ywkd2e4ko5x3/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o -MF CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o.d -o CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o -c /cluster/home/tamaurer/lehrfempp/lehrfempp/projects/dgfem/fe_experiments.cc

projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.i"
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem && /cluster/apps/gcc-4.8.5/gcc-9.3.0-mjmsvu3362usxjr7tas4ywkd2e4ko5x3/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cluster/home/tamaurer/lehrfempp/lehrfempp/projects/dgfem/fe_experiments.cc > CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.i

projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.s"
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem && /cluster/apps/gcc-4.8.5/gcc-9.3.0-mjmsvu3362usxjr7tas4ywkd2e4ko5x3/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cluster/home/tamaurer/lehrfempp/lehrfempp/projects/dgfem/fe_experiments.cc -o CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.s

# Object files for target projects.dgfem.meshes
projects_dgfem_meshes_OBJECTS = \
"CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o"

# External object files for target projects.dgfem.meshes
projects_dgfem_meshes_EXTERNAL_OBJECTS =

projects/dgfem/projects.dgfem.meshes: projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/fe_experiments.cc.o
projects/dgfem/projects.dgfem.meshes: projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/build.make
projects/dgfem/projects.dgfem.meshes: /cluster/home/tamaurer/.hunter/_Base/093bd9a/872b3a8/84133a3/Install/lib/libboost_program_options-mt-d-x64.a
projects/dgfem/projects.dgfem.meshes: lib/lf/assemble/liblf.assemble.a
projects/dgfem/projects.dgfem.meshes: lib/lf/base/liblf.base.a
projects/dgfem/projects.dgfem.meshes: lib/lf/geometry/liblf.geometry.a
projects/dgfem/projects.dgfem.meshes: lib/lf/io/liblf.io.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/liblf.mesh.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/hybrid2d/liblf.mesh.hybrid2d.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/polytopic2d/liblf.mesh.polytopic2d.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/test_utils/liblf.mesh.test_utils.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/utils/liblf.mesh.utils.a
projects/dgfem/projects.dgfem.meshes: lib/lf/refinement/liblf.refinement.a
projects/dgfem/projects.dgfem.meshes: lib/lf/fe/liblf.fe.a
projects/dgfem/projects.dgfem.meshes: lib/lf/uscalfe/liblf.uscalfe.a
projects/dgfem/projects.dgfem.meshes: /cluster/home/tamaurer/.hunter/_Base/093bd9a/872b3a8/84133a3/Install/lib64/libgtest_maind.a
projects/dgfem/projects.dgfem.meshes: /cluster/home/tamaurer/.hunter/_Base/093bd9a/872b3a8/84133a3/Install/lib64/libgtestd.a
projects/dgfem/projects.dgfem.meshes: lib/lf/io/liblf.io.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/polytopic2d/liblf.mesh.polytopic2d.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/hybrid2d/liblf.mesh.hybrid2d.a
projects/dgfem/projects.dgfem.meshes: lib/lf/fe/liblf.fe.a
projects/dgfem/projects.dgfem.meshes: lib/lf/assemble/liblf.assemble.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/liblf.mesh.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/utils/liblf.mesh.utils.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/liblf.mesh.a
projects/dgfem/projects.dgfem.meshes: lib/lf/mesh/utils/liblf.mesh.utils.a
projects/dgfem/projects.dgfem.meshes: lib/lf/geometry/liblf.geometry.a
projects/dgfem/projects.dgfem.meshes: lib/lf/quad/liblf.quad.a
projects/dgfem/projects.dgfem.meshes: lib/lf/base/liblf.base.a
projects/dgfem/projects.dgfem.meshes: /cluster/home/tamaurer/.hunter/_Base/093bd9a/872b3a8/84133a3/Install/lib64/libspdlogd.a
projects/dgfem/projects.dgfem.meshes: /cluster/home/tamaurer/.hunter/_Base/093bd9a/872b3a8/84133a3/Install/lib64/libfmtd.a
projects/dgfem/projects.dgfem.meshes: projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cluster/home/tamaurer/lehrfempp/lehrfempp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable projects.dgfem.meshes"
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/projects.dgfem.meshes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/build: projects/dgfem/projects.dgfem.meshes
.PHONY : projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/build

projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/clean:
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem && $(CMAKE_COMMAND) -P CMakeFiles/projects.dgfem.meshes.dir/cmake_clean.cmake
.PHONY : projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/clean

projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/depend:
	cd /cluster/home/tamaurer/lehrfempp/lehrfempp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cluster/home/tamaurer/lehrfempp/lehrfempp /cluster/home/tamaurer/lehrfempp/lehrfempp/projects/dgfem /cluster/home/tamaurer/lehrfempp/lehrfempp/build /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem /cluster/home/tamaurer/lehrfempp/lehrfempp/build/projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : projects/dgfem/CMakeFiles/projects.dgfem.meshes.dir/depend

