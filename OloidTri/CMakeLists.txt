
cmake_minimum_required(VERSION 2.8)

PROJECT(OloidTri)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(OloidTri MACOSX_BUNDLE OloidTri)

if(VTK_LIBRARIES)
  target_link_libraries(OloidTri ${VTK_LIBRARIES})
else()
  target_link_libraries(OloidTri vtkHybrid vtkWidgets)
endif()
