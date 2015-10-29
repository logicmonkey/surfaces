
cmake_minimum_required(VERSION 2.8)

PROJECT(FoloidRuled)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FoloidRuled MACOSX_BUNDLE FoloidRuled)

if(VTK_LIBRARIES)
  target_link_libraries(FoloidRuled ${VTK_LIBRARIES})
else()
  target_link_libraries(FoloidRuled vtkHybrid vtkWidgets)
endif()
