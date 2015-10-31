
cmake_minimum_required(VERSION 2.8)

PROJECT(Oloid)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

add_executable(Oloid MACOSX_BUNDLE Oloid)

if(VTK_LIBRARIES)
  target_link_libraries(Oloid ${VTK_LIBRARIES})
else()
  target_link_libraries(Oloid vtkHybrid vtkWidgets)
endif()
