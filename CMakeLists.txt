#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DBDecay)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIR}")
execute_process(COMMAND root-config --incdir OUTPUT_VARIABLE ROOT_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
#set(CRY_INCLUDE_DIR "/mnt/hgfs/ShareE/GEANT4/cry_v1.7/cry_v1.7")

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(DBDecay DBDecay.cc ${sources} ${headers})
target_link_libraries(DBDecay ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Link the CRY package
#
#set(CRY_STATIC_LIB "/mnt/hgfs/ShareE/GEANT4/cry_v1.7/cry_v1.7/lib/libCRY.a")
#target_link_libraries(DBDecay ${CRY_STATIC_LIB})

#----------------------------------------------------------------------------
# Copy all scripts from macFiles folder to the build directory, i.e. the directory in which we
# build DBDecay. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB DBDecay_SCRIPTS "${PROJECT_SOURCE_DIR}/macFiles/*.mac")

foreach(_script ${DBDecay_SCRIPTS})
  get_filename_component(_script_name ${_script} NAME)

  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/${_script_name}
    COPYONLY
  )
endforeach()

#configure_file(
#  ${PROJECT_SOURCE_DIR}/macFiles/mergeResults.cc
#  ${PROJECT_BINARY_DIR}/resultFiles/mergeResults.cc
#  COPYONLY
#)

#----------------------------------------------------------------------------
# create the result folder
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/resultFiles)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/printfile)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DBDecay DESTINATION bin)

