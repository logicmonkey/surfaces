
cmake_minimum_required(VERSION 2.8)

PROJECT(FoloidMesh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(-DDEBUG)

add_executable(FoloidMesh MACOSX_BUNDLE FoloidMesh)

if(VTK_LIBRARIES)
  target_link_libraries(FoloidMesh ${VTK_LIBRARIES})
else()
  target_link_libraries(FoloidMesh vtkHybrid vtkWidgets)
endif()
