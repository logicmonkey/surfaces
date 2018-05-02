Package requirements for a build on Arch Linux with VTK version 8.x.x
sudo pacman -Syu
sudo pacman -S vtk cmake qt5-base qt5-x11extras glew hdf5 netcdf-cxx proj openmpi gdal unixodbc

Example:
cd rOloid
mkdir build
cd build
cmake .. && make
