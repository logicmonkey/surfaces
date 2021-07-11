
cmake_minimum_required(VERSION 2.8)

PROJECT(TrefoilKnot)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

set(SRC
  TrefoilKnot.cxx
  main.cxx
  )

set(INC
  TrefoilKnot.h
  )

add_executable(TrefoilKnot MACOSX_BUNDLE ${SRC} ${INC})

if(VTK_LIBRARIES)
  target_link_libraries(TrefoilKnot ${VTK_LIBRARIES})
else()
  target_link_libraries(TrefoilKnot vtkHybrid vtkWidgets)
endif()
