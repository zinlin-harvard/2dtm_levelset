cgal_create_CMakeLists -s convert_gmsh_to_xml
cmake -DCGAL_DIR=$HOME/MyLocal/CGAL-4.7-Installed .
make
rm -rf Makefile cmake_install.cmake CMake*
