
cmake_minimum_required(VERSION 2.8)

PROJECT(FauxloidMesh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

add_executable(FauxloidMesh MACOSX_BUNDLE FauxloidMesh)

if(VTK_LIBRARIES)
  target_link_libraries(FauxloidMesh ${VTK_LIBRARIES})
else()
  target_link_libraries(FauxloidMesh vtkHybrid vtkWidgets)
endif()
