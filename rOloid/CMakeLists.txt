
cmake_minimum_required(VERSION 2.8)

PROJECT(rOloid)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

add_executable(rOloid MACOSX_BUNDLE rOloid)

if(VTK_LIBRARIES)
  target_link_libraries(rOloid ${VTK_LIBRARIES})
else()
  target_link_libraries(rOloid vtkHybrid vtkWidgets)
endif()
