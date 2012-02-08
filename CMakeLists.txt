#******************************************************************************
#
# MantaFlow fluid solver framework
#
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL) 
# http://www.gnu.org/licenses
#
#******************************************************************************

project (MantaFlow)
cmake_minimum_required(VERSION 2.8) # 2.8 needed for CUDA

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/source/cmake/")
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
SET(VERBOSE 1) 

#******************************************************************************
# setup vars

IF(NOT CMAKE_BUILD_TYPE) 
	SET(CMAKE_BUILD_TYPE "Release") 
ELSE()
	MESSAGE("Build-type: '${CMAKE_BUILD_TYPE}'")
ENDIF()

# compilation versions
OPTION(MT "Use multi-threaded kernels (TBB)" OFF)
OPTION(DEBUG "Enable debug compilation" OFF)
OPTION(USE64 "Switch 64bit support" ON)
OPTION(GUI "Compile with GUI (requires QT)" OFF)
OPTION(PREPDEBUG "Debug generated files" OFF) # This will beautify generated files, and link to them for compiler errors instead of the original sources
OPTION(CUDA "Compile with CUDA plugins" OFF)
OPTION(MESHCODE "Compile with meshcode" OFF)

# make sure debug settings match...
IF(NOT DEBUG)
	IF(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
		SET(DEBUG 1) 
	ENDIF()
ENDIF()
IF(DEBUG)
	SET(CMAKE_BUILD_TYPE "Debug")
    add_definitions ( -DDEBUG )
ENDIF()

MESSAGE("")
MESSAGE("Options - "
    " -DMT='${MT}' "
    " -DDEBUG='${DEBUG}' "
    " -DPREPDEBUG='${PREPDEBUG}' "
    " -DUSE64='${USE64}' "
    " -DGUI='${GUI}' "
    " -DCUDA='${CUDA}' "
    )
MESSAGE("")

#******************************************************************************
# Pre-processor

# compile prep
SET(SOURCES 
    source/preprocessor/main.cpp
    source/preprocessor/parse.cpp
    source/preprocessor/merge.cpp
    source/preprocessor/replace.cpp
    source/preprocessor/tokenize.cpp
)
add_executable(prep ${SOURCES})
if (NOT WIN32)
	set_target_properties(prep PROPERTIES COMPILE_FLAGS "-Wall -O2")
endif()

#******************************************************************************
# Setup main project

SET(F_LIBS "" )
SET(F_LIB_PATHS)
SET(F_LINKADD "") # additional linker flags, not a list
set(PP_PATH "pp")

# need pre-processing
SET(PP_SOURCES
    source/fluidsolver.cpp
    source/conjugategrad.cpp
    source/grid.cpp
    source/particle.cpp
    source/levelset.cpp
    source/fastmarch.cpp
    source/shapes.cpp
    source/mesh.cpp
    source/general.cpp
    source/fileio.cpp   
    source/noisefield.cpp
    source/kernel.cpp
    source/plugin/pressure.cpp
    source/plugin/advection.cpp
    source/plugin/extforces.cpp
    source/python/defines.py
)

SET(PP_HEADERS
    source/commonkernels.h
    source/conjugategrad.h
    source/fastmarch.h
    source/fluidsolver.h
    source/grid.h
    source/mesh.h
    source/particle.h
    source/levelset.h
    source/shapes.h
    source/noisefield.h
)

# no pre-processing needed
set(SOURCES
    source/pwrapper/pymain.cpp
    source/pwrapper/pclass.cpp
    source/pwrapper/pvec3.cpp
    source/pwrapper/pconvert.cpp
    source/util/vectorbase.cpp
)
 
SET(HEADERS
    source/pwrapper/pclass.h
    source/pwrapper/pconvert.h
    source/util/vectorbase.h
    source/general.h
    source/kernel.h
    source/fileio.h
)

if (GUI)
    # need QT preprocessor
    set(QT_HEADERS 
        source/gui/mainwindow.h
        source/gui/glwidget.h
        source/gui/painter.h
        source/gui/meshpainter.h
        source/gui/qtmain.h 
        source/gui/customctrl.h
        source/gui/particlepainter.h
    )
    list(APPEND PP_SOURCES
        source/gui/customctrl.cpp
        source/gui/mainwindow.cpp
        source/gui/glwidget.cpp
        source/gui/customctrl.cpp
        source/gui/painter.cpp
        source/gui/meshpainter.cpp
        source/gui/particlepainter.cpp
        source/gui/qtmain.cpp
    )
endif()

# CUDA sources
if (CUDA)
    list(APPEND PP_SOURCES 
        source/cuda/meshtools.cu
        source/cuda/particle.cu
    )
endif()

if(MESHCODE)
    add_definitions(-DMESHCODE)
	list(APPEND PP_HEADERS
		source/edgecollapse.h
		source/vortexsheet.h
		source/vortexpart.h
    )
    list(APPEND PP_SOURCES
        source/test.cpp    
        source/edgecollapse.cpp
        source/plugin/vortexplugins.cpp
		source/vortexsheet.cpp    
		source/vortexpart.cpp    
        source/plugin/meshplugins.cpp    
    )
    if (CUDA)
        list(APPEND PP_SOURCES 
            source/cuda/buoyancy.cu
            source/cuda/turbulence.cu
            source/cuda/curlnoise.cu        
            source/plugin/turbulence.cpp        
        )
    endif()
endif()

# include dirs
SET(INCLUDE_PATHS
   source/pwrapper
   source/util
   source/cuda
   ${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source
   source
)

# 64 bit config
if (USE64)
    SET(F_LINKADD "${F_LINKADD} -m64 ")
    ADD_DEFINITIONS( -m64 )
endif()

# Multithreading
if(MT)
    add_definitions( -DMT)
    if (DEBUG)
        add_definitions( -DTBB_USE_DEBUG=1 )
        list(APPEND F_LIBS tbb)
    else()
        list(APPEND F_LIBS tbb)
    endif()
	if (WIN32)
		find_package(TBB REQUIRED)
	endif()
endif()

#******************************************************************************
# add a target to generate API documentation with Doxygen

find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DX_PATH "doxy")
    foreach(it ${PP_SOURCES} ${PP_HEADERS} ${SOURCES} ${QT_HEADERS} ${PP_HEADERS})
        get_filename_component(CURPATH ${it} PATH)
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${DX_PATH}/${CURPATH}")
        set(CURDX "${DX_PATH}/${it}")
        string(REPLACE "/" "_" TGT ${CURDX})
        add_custom_command(OUTPUT ${TGT} 
                        COMMAND ./prep docgen "${CMAKE_CURRENT_SOURCE_DIR}/${it}" "${CURDX}"
                        DEPENDS prep
                        IMPLICIT_DEPENDS CXX ${it}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        list(APPEND TGLIST ${TGT})
    endforeach(it)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc 
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${TGLIST}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

#******************************************************************************
# Link libraries

# Python 3.2
if(APPLE)
    #hack: cmake/mac can't find python3 for some reason
    list(APPEND INCLUDE_PATHS /Library/Frameworks/Python.framework/Versions/3.2/include/python3.2m)
    list(APPEND F_LIB_PATHS /Library/Frameworks/Python.framework/Versions/3.2/lib)
    list(APPEND F_LIBS "python3.2")
else()
    set(Python_ADDITIONAL_VERSIONS 3.1 3.2)
    find_package(PythonLibs REQUIRED)
    list(APPEND INCLUDE_PATHS ${PYTHON_INCLUDE_DIRS})
    list(APPEND F_LIBS ${PYTHON_LIBRARIES})
endif()

# Z compression
include(FindZLIB)
LIST(APPEND F_LIBS "z")
LIST(APPEND INCLUDE_PATHS ${ZLIB_INCLUDE_DIR})

# CUDA
if(CUDA)
    add_definitions( -DCUDA=1 )
    find_package(CUDA QUIET REQUIRED)
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
    set(CUDA_64_BIT_DEVICE_CODE ON)
    
    if(DEBUG)
        set(CUDA_NVCC_FLAGS -g;-ccbin;gcc-4.5;-DDEBUG;-keep;-arch=sm_20;--maxrregcount=31;--compiler-options;-Wall)
    else()
        set(CUDA_NVCC_FLAGS -ccbin;gcc-4.5;--use_fast_math;-arch=sm_20;-DNDEBUG;--compiler-options;-Wall;-O3)       
    endif()    
endif()

#******************************************************************************
# apply preprocessor

set(PP_REGISTRY "${PP_PATH}/ppreg.cpp")
set(PP_REGS)
set(PP_PREPD "0")
if (PREPDEBUG)
    set(PP_PREPD "1")
endif()
set(PP_USEMT "0")
if (MT)
    set(PP_USEMT "1")
endif()
foreach(it ${PP_SOURCES} ${PP_HEADERS} ${QT_HEADERS})
    get_filename_component(CURPATH ${it} PATH)
    #get_filename_component(CUREXT ${it} EXT)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${CURPATH}")
    set(CURPP "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${it}")
    # preprocessor
    add_custom_command(OUTPUT ${CURPP}
                       COMMAND ./prep generate ${PP_PREPD} ${PP_USEMT} "${CMAKE_CURRENT_SOURCE_DIR}/${it}" "${CURPP}"
                       DEPENDS prep
                       IMPLICIT_DEPENDS CXX ${it}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    list(APPEND PP_REGS "${CURPP}.reg")
    list(APPEND SOURCES ${CURPP})    
endforeach(it)

# merge reg files
add_custom_command(OUTPUT ${PP_REGISTRY}
                  COMMAND ./prep merge ${PP_REGISTRY} ${PP_REGS}
                  DEPENDS prep ${SOURCES}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_source_files_properties(${PP_REGISTRY} PROPERTIES GENERATED 1)

#******************************************************************************
# QT for GUI

if(GUI)
    find_package(Qt4 REQUIRED)
    set(QT_USE_QTOPENGL TRUE)

    # remap
    set(QT_REMAP)
    foreach(it ${QT_HEADERS})
        list(APPEND QT_REMAP "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${it}")
    endforeach(it)
    
    qt4_wrap_cpp(MOC_OUTFILES ${QT_REMAP} )
    qt4_add_resources(QT_RES resources/res.qrc )
    include(${QT_USE_FILE})
    add_definitions(-DGUI=1)
    add_definitions(${QT_DEFINITIONS})
    list(APPEND INCLUDE_PATHS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source/gui source/gui)
    list(APPEND F_LIBS ${QT_LIBRARIES})
    list(APPEND SOURCES ${MOC_OUTFILES} ${QT_RES})
    
    if (APPLE)
        SET(F_LINKADD "${F_LINKADD} -framework OpenGL ")
    else()
        list(APPEND F_LIBS "GL")
    endif()
endif()

#******************************************************************************
# setup executable

# compiler flags
IF(NOT WIN32)
    IF(DEBUG)
        add_definitions( -DDEBUG )
        add_definitions( -O0 -fno-inline )
        # additional c++ warning
        SET(CMAKE_CXX_FLAGS_DEBUG " -ggdb  -Wall -Wno-reorder")
        
        # stricter: no optimizations and inlining, comment out if needed...
        #SET(CMAKE_CXX_FLAGS_DEBUG " -O0 -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}   ")              
    ELSE()
        # non-debug, optimized version
        add_definitions( -O3 )
    ENDIF()
ELSE()
    # get rid of some MSVC warnings
    add_definitions( /wd4018 /wd4146 /wd4800 )
ENDIF()

SET(EXECCMD manta)

SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${F_LINKADD} ")
SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${F_LINKADD} ")

include_directories( ${INCLUDE_PATHS}) 
link_directories( ${F_LIB_PATHS} ) 

if(CUDA)
    #cuda_include_directories(pp/source/cuda)
    cuda_add_executable( ${EXECCMD} ${SOURCES} ${PP_REGISTRY})
    target_link_libraries( ${EXECCMD} ${F_LIBS} )
    cuda_build_clean_target()
else()
    add_executable(${EXECCMD} ${SOURCES} ${PP_REGISTRY})
    target_link_libraries( ${EXECCMD} ${F_LIBS} )
endif()

# debug
MESSAGE("DEBUG Flag-Summary - Includes: '${INCLUDE_PATHS}' | Libs: '${F_LIBS}' | LibPaths: '${F_LIB_PATHS}' ")


