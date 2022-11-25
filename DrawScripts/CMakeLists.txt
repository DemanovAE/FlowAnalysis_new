# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(FemtoDstAnalyzer)

if (ROOT_CMAKE)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else (ROOT_CMAKE)
  set(ROOT_PREINSTALLED "/usr/lib64/Geant4-10.3.0/Modules")
  list(APPEND CMAKE_MODULE_PATH ${ROOT_PREINSTALLED})
endif (ROOT_CMAKE)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net MathMore EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Define _VANILLA_ROOT_ variable needed for the project
#add_definitions(-D_VANILLA_ROOT_)

add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set(CMAKE_CXX_STANDARD 17)


set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
)

set(DrawPicture_LIBRARY_h_files
  ${CMAKE_SOURCE_DIR}/MyDrawObject.h
  ${CMAKE_SOURCE_DIR}/Draw_Picture_new.h
  ${CMAKE_SOURCE_DIR}/AnalysisFunction.h
)

set(DrawPicture_LIBRARY_cxx_files
  ${CMAKE_SOURCE_DIR}/MyDrawObject.cxx
  ${CMAKE_SOURCE_DIR}/Draw_Picture_new.cxx
  ${CMAKE_SOURCE_DIR}/AnalysisFunction.cxx
)

set(DrawPicture_LinkDef
  ${CMAKE_SOURCE_DIR}/DrawPict.LinkDef.h
)

include_directories(${INCLUDE_DIRECTORIES})

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G__DrawPicture
  ${DrawPicture_LIBRARY_h_files}
  LINKDEF ${DrawPicture_LinkDef}
)

#---Compile library
add_library(DrawPicture SHARED ${DrawPicture_LIBRARY_cxx_files} G__DrawPicture.cxx)
target_link_libraries(DrawPicture ${ROOT_LIBRARIES})

set(DrawPicture_INCLUDE_LIBRARIES
  ${ROOT_LIBRARIES}
  DrawPicture
)

#---Make executable

#add_executable(Draw DrawV2.cpp)
#target_link_libraries(Draw ${DrawPicture_INCLUDE_LIBRARIES})