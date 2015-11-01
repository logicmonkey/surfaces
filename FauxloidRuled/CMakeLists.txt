
cmake_minimum_required(VERSION 2.8)

PROJECT(FauxloidRuled)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FauxloidRuled MACOSX_BUNDLE FauxloidRuled)

if(VTK_LIBRARIES)
  target_link_libraries(FauxloidRuled ${VTK_LIBRARIES})
else()
  target_link_libraries(FauxloidRuled vtkHybrid vtkWidgets)
endif()
