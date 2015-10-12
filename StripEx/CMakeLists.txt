
cmake_minimum_required(VERSION 2.8)

PROJECT(StripEx)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(StripEx MACOSX_BUNDLE StripEx)

if(VTK_LIBRARIES)
  target_link_libraries(StripEx ${VTK_LIBRARIES})
else()
  target_link_libraries(StripEx vtkHybrid vtkWidgets)
endif()
