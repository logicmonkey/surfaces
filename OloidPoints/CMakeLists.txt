
cmake_minimum_required(VERSION 2.8)

PROJECT(OloidPoints)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

add_executable(OloidPoints MACOSX_BUNDLE OloidPoints)

if(VTK_LIBRARIES)
  target_link_libraries(OloidPoints ${VTK_LIBRARIES})
else()
  target_link_libraries(OloidPoints vtkHybrid vtkWidgets)
endif()
