
cmake_minimum_required(VERSION 2.8)

PROJECT(OloidMesh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(OloidMesh MACOSX_BUNDLE OloidMesh)

if(VTK_LIBRARIES)
  target_link_libraries(OloidMesh ${VTK_LIBRARIES})
else()
  target_link_libraries(OloidMesh vtkHybrid vtkWidgets)
endif()
