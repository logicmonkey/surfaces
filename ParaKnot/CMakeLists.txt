
cmake_minimum_required(VERSION 2.8)

PROJECT(ParaKnot)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

set(SRC
  ParaKnot.cxx
  main.cxx
  )

set(INC
  ParaKnot.h
  )

add_executable(ParaKnot MACOSX_BUNDLE ${SRC} ${INC})

if(VTK_LIBRARIES)
  target_link_libraries(ParaKnot ${VTK_LIBRARIES})
else()
  target_link_libraries(ParaKnot vtkHybrid vtkWidgets)
endif()
