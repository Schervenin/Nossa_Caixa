#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(NossaCaixa1.0)

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
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(NossaCaixa1.0 NossaCaixa1.0.cc ${sources} ${headers})
target_link_libraries(NossaCaixa1.0 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

set(Calor5_SCRIPTS
    initInter.mac
#    initInterAll.mac 
    Calor_vis.mac 
#    Guide_NossaCaixa1.0.pdf
)

# Directories - uncomment for copy in the build directory
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples)

#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Gamma)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Gamma/AngularDistribution)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Gamma/Isotropic)

#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/IntrinsicRadiation)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/IntrinsicRadiation/LaBr3)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/IntrinsicRadiation/LYSO)

#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Isotope)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Isotope/Co60)
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/Isotope/DopplerShift)

#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Examples/PhaseSpace)

file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/DetectorsPositionFiles)

file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Wkq_files)

# Files
#file(GLOB _ex_gamma_ang ${PROJECT_SOURCE_DIR}/Examples/Gamma/AngularDistribution/*)
#file(GLOB _ex_gamma_iso ${PROJECT_SOURCE_DIR}/Examples/Gamma/Isotropic/*)

#file(GLOB _ex_intrinsic_LaBr ${PROJECT_SOURCE_DIR}/Examples/IntrinsicRadiation/LaBr3/*)
#file(GLOB _ex_intrinsic_LYSO ${PROJECT_SOURCE_DIR}/Examples/IntrinsicRadiation/LYSO/*)

#file(GLOB _ex_isot_Co60 ${PROJECT_SOURCE_DIR}/Examples/Isotope/Co60/*)
#file(GLOB _ex_isot_doppler ${PROJECT_SOURCE_DIR}/Examples/Isotope/DopplerShift/*)

#file(GLOB _ex_PhaseSpace ${PROJECT_SOURCE_DIR}/Examples/PhaseSpace/*)

#file(GLOB _ex_DetPos ${PROJECT_SOURCE_DIR}/DetectorsPositionFiles/*)

#file(GLOB _ex_Wkq ${PROJECT_SOURCE_DIR}/Wkq_files/*)

# Copy
foreach(_script ${Calor5_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#foreach(_ex_gamma_angfile ${_ex_gamma_ang})
#   configure_file( ${_ex_gamma_angfile} ${PROJECT_BINARY_DIR}/Examples/Gamma/AngularDistribution COPYONLY)    
#endforeach()

#foreach(_ex_gamma_isofile ${_ex_gamma_iso})
#   configure_file( ${_ex_gamma_isofile} ${PROJECT_BINARY_DIR}/Examples/Gamma/Isotropic COPYONLY)    
#endforeach()

#foreach(_ex_intrinsic_LaBrfile ${_ex_intrinsic_LaBr})
#   configure_file( ${_ex_intrinsic_LaBrfile} ${PROJECT_BINARY_DIR}/Examples/IntrinsicRadiation/LaBr3 COPYONLY)    
#endforeach()

#foreach(_ex_intrinsic_LYSOfile ${_ex_intrinsic_LYSO})
#   configure_file( ${_ex_intrinsic_LYSOfile} ${PROJECT_BINARY_DIR}/Examples/IntrinsicRadiation/LYSO COPYONLY)    
#endforeach()

#foreach(_ex_isot_Co60file ${_ex_isot_Co60})
#   configure_file( ${_ex_isot_Co60file} ${PROJECT_BINARY_DIR}/Examples/Isotope/Co60 COPYONLY)    
#endforeach()

#foreach(_ex_isot_dopplerfile ${_ex_isot_doppler})
#   configure_file( ${_ex_isot_dopplerfile} ${PROJECT_BINARY_DIR}/Examples/Isotope/DopplerShift COPYONLY)    
#endforeach()

#foreach(_ex_PhaseSpacefile ${_ex_PhaseSpace})
#   configure_file( ${_ex_PhaseSpacefile} ${PROJECT_BINARY_DIR}/Examples/PhaseSpace COPYONLY)    
#endforeach()


foreach(_ex_DetPosfile ${_ex_DetPos})
   configure_file( ${_ex_DetPosfile} ${PROJECT_BINARY_DIR}/DetectorsPositionFiles COPYONLY)    
endforeach()

foreach(_ex_Wkqfile ${_ex_Wkq})
   configure_file( ${_ex_Wkqfile} ${PROJECT_BINARY_DIR}/Wkq_files COPYONLY)    
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(NossaCaixa1.0 DEPENDS exampleB1)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS NossaCaixa1.0 DESTINATION bin)


