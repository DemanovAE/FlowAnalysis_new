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

set(ST_FEMTO_DST_LIB libStFemtoDst.so)
set(ST_FEMTO_DST_CUSTOM_DIR $ENV{ST_FEMTO_DST_BUILD_DIR})
set(ST_FEMTO_DST_CUSTOM_INC_DIR $ENV{ST_FEMTO_DST_INC_DIR})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set(CMAKE_CXX_STANDARD 17)

find_library(ST_FEMTO_DST REQUIRED
	NAMES ${ST_FEMTO_DST_LIB}
	HINTS ${ST_FEMTO_DST_CUSTOM_DIR})

set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}
  ${ST_FEMTO_DST_CUSTOM_INC_DIR}
  ${ROOT_INCLUDE_DIRS}
)

set(FlowAnalysis_LIBRARY_h_files
  ${CMAKE_SOURCE_DIR}/FlowEtaSubEP.h
  ${CMAKE_SOURCE_DIR}/FlowScalarProduct.h
  ${CMAKE_SOURCE_DIR}/QVector.h
  ${CMAKE_SOURCE_DIR}/PIDflow.h
  ${CMAKE_SOURCE_DIR}/StEpdGeom.h
  ${CMAKE_SOURCE_DIR}/QAFemtoDst2.h
)

set(FlowAnalysis_LIBRARY_cxx_files
  ${CMAKE_SOURCE_DIR}/FlowEtaSubEP.cxx
  ${CMAKE_SOURCE_DIR}/FlowScalarProduct.cxx
  ${CMAKE_SOURCE_DIR}/QVector.cxx
  ${CMAKE_SOURCE_DIR}/PIDflow.cxx
  ${CMAKE_SOURCE_DIR}/StEpdGeom.cxx
  ${CMAKE_SOURCE_DIR}/QAFemtoDst2.cxx
)

set(FlowAnalysis_LinkDef
  ${CMAKE_SOURCE_DIR}/FlowAnalysis.LinkDef.h
)

include_directories(${INCLUDE_DIRECTORIES})

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G__FlowAnalysis
  ${FlowAnalysis_LIBRARY_h_files}
  LINKDEF ${FlowAnalysis_LinkDef}
)

#---Compile library
add_library(FlowAnalysis SHARED ${FlowAnalysis_LIBRARY_cxx_files} G__FlowAnalysis.cxx)
target_link_libraries(FlowAnalysis ${ROOT_LIBRARIES})

set(FemtoDstAnalyzer_INCLUDE_LIBRARIES
  ${ROOT_LIBRARIES}
  ${ST_FEMTO_DST}
  FlowAnalysis
)

#---Make executable

add_executable(FemtoDstAnalyzer_New main_Hadrons.cpp)
target_link_libraries(FemtoDstAnalyzer_New ${FemtoDstAnalyzer_INCLUDE_LIBRARIES})
add_executable(FemtoDstAnalyzer_QA main_QA.cpp)
target_link_libraries(FemtoDstAnalyzer_QA ${FemtoDstAnalyzer_INCLUDE_LIBRARIES})
add_executable(combPIDFit main_FitCombPID.cpp)
target_link_libraries(combPIDFit ${FemtoDstAnalyzer_INCLUDE_LIBRARIES})