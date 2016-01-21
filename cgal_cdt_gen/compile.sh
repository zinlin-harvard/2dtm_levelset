cgal_create_CMakeLists -s gen_optimized_mesh
cmake -DCGAL_DIR=$HOME/MyLocal/CGAL-4.7-Installed .
make
rm -rf Makefile cmake_install.cmake CMake*
