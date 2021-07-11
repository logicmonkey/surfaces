
cmake_minimum_required(VERSION 2.8)

PROJECT(Figure8Torus)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

set(SRC
  Figure8Torus.cxx
  main.cxx
  )

set(INC
  Figure8Torus.h
  )

add_executable(Figure8Torus MACOSX_BUNDLE ${SRC} ${INC})

if(VTK_LIBRARIES)
  target_link_libraries(Figure8Torus ${VTK_LIBRARIES})
else()
  target_link_libraries(Figure8Torus vtkHybrid vtkWidgets)
endif()
