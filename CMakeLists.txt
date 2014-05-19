#******************************************************************************
#
# MantaFlow fluid solver framework
#
# Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
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

SET(MANTAVERSION "0.6")

#******************************************************************************
# Default paths
# - CMake's path finder is completely useless for Qt5 + Python on Win64

IF(WIN32)
    set(WIN_QT_PATH "C:/Qt/Qt5.2.1/5.2.1/msvc2012_64_opengl") # qt5/win64
    # set(WIN_QT_PATH "C:/Qt/4.8.6") # qt4/win32
    set(WIN_PYTHON_PATH "C:/Python34")
	set(WIN_PYTHON_VERSION "34")
	set(CMAKE_LIBRARY_PATH "C:/Program Files (x86)/Windows Kits/8.0/Lib/win8/um/x64")
    SET(CMAKE_PREFIX_PATH ${WIN_QT_PATH})
ENDIF()

IF(APPLE)
    SET(CMAKE_PREFIX_PATH "/usr/local/Cellar/qt5/5.2.1/") #  mac/homebrew
ENDIF()

#******************************************************************************
# setup default params

IF (CMAKE_SIZEOF_VOID_P EQUAL 8) 
	SET(64BIT TRUE)
ELSE()
	SET(64BIT FALSE)
ENDIF()

IF(NOT CMAKE_BUILD_TYPE) 
    SET(CMAKE_BUILD_TYPE "Release") 
ELSE()
    MESSAGE("Build-type: '${CMAKE_BUILD_TYPE}'")
ENDIF()

# compilation versions
OPTION(DEBUG "Enable debug compilation" OFF)
OPTION(GUI "Compile with GUI (requires QT)" OFF)
OPTION(TBB "Use multi-thread kernels using Intels TBB" OFF)
OPTION(OPENMP "Use multi-thread kernels using OpenMP" OFF)
OPTION(CUDA "Compile with CUDA plugins" OFF)
OPTION(PREPDEBUG "Debug generated files" OFF) # This will beautify generated files, and link to them for compiler errors instead of the original sources
OPTION(DOUBLEPRECISION "Compile with double floating point precision" OFF)

#check consistency of MT options
set(MT OFF)
set(MT_TYPE "NONE")
if (TBB)
    set (MT_TYPE "TBB")
    set (MT ON)
endif()
if (OPENMP)
    set (MT_TYPE "OPENMP")
    set (MT ON)
endif()
if (TBB AND OPENMP)
    message(FATAL_ERROR "Cannot activate both OPENMP and TBB")
endif()

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

MESSAGE(STATUS "")
MESSAGE(STATUS "Options - "
    " -DDEBUG='${DEBUG}' "
    " -DGUI='${GUI}' "
    " -DTBB='${TBB}' "
    " -DOPENMP='${OPENMP}' "
    " -DPREPDEBUG='${PREPDEBUG}' "
    " -DCUDA='${CUDA}' "
    " -DDOUBLEPRECISION='${DOUBLEPRECISION}' "
    )
MESSAGE(STATUS "Multithreading type : ${MT_TYPE}")
MESSAGE(STATUS "")

#******************************************************************************
# Pre-processor

# compile prep
SET(SOURCES 
    source/preprocessor/main.cpp
    source/preprocessor/code.cpp
    source/preprocessor/tokenize.cpp
    source/preprocessor/parse.cpp
    source/preprocessor/util.cpp
    source/preprocessor/merge.cpp
    source/preprocessor/codegen_python.cpp
    source/preprocessor/codegen_kernel.cpp    
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
SET(SILENT_SOURCES)

# need pre-processing
SET(PP_SOURCES
    source/fluidsolver.cpp
    source/conjugategrad.cpp
    source/grid.cpp
    source/levelset.cpp
    source/fastmarch.cpp
    source/shapes.cpp
    source/mesh.cpp
    source/general.cpp
    source/particle.cpp
    source/flip.cpp
    source/movingobs.cpp
    source/fileio.cpp   
    source/noisefield.cpp
    source/kernel.cpp
    source/vortexsheet.cpp
    #source/vortexfilament.cpp
    source/vortexpart.cpp
    source/turbulencepart.cpp
    source/edgecollapse.cpp
    source/plugin/pressure.cpp
    source/plugin/initplugins.cpp
    source/plugin/advection.cpp
    source/plugin/extforces.cpp
    source/plugin/kepsilon.cpp
    source/plugin/meshplugins.cpp
    source/plugin/vortexplugins.cpp
    source/plugin/waveletturbulence.cpp
    source/python/defines.py
    source/test.cpp
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
    source/vortexsheet.h
    #source/vortexfilament.h
    source/general.h
    source/kernel.h
    source/flip.h
    source/movingobs.h
    source/fileio.h
    source/edgecollapse.h
    source/vortexpart.h
    source/turbulencepart.h
    #source/test.h    
)

# no pre-processing needed
set(NOPP_SOURCES
    source/pwrapper/pymain.cpp
    source/pwrapper/pclass.cpp
    source/pwrapper/pvec3.cpp
    source/pwrapper/pconvert.cpp
    source/pwrapper/registry.cpp
    source/util/vectorbase.cpp
)
 
SET(NOPP_HEADERS
    source/pwrapper/pythonInclude.h
    source/pwrapper/pclass.h
    source/pwrapper/registry.h
    source/pwrapper/pconvert.h
    source/util/integrator.h
    source/util/vectorbase.h
    source/util/quaternion.h
    source/util/interpol.h
    source/util/mcubes.h
    source/util/randomstream.h
    source/util/solvana.h    
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
    set(QT_SOURCES
        source/gui/customctrl.cpp
        source/gui/mainwindow.cpp
        source/gui/glwidget.cpp
        source/gui/customctrl.cpp
        source/gui/painter.cpp
        source/gui/meshpainter.cpp
        source/gui/particlepainter.cpp
        source/gui/qtmain.cpp
    )
	list(APPEND PP_SOURCES ${QT_SOURCES})
	list(APPEND PP_HEADERS ${QT_HEADERS})
endif()

# CUDA sources
if (CUDA)
    list(APPEND PP_SOURCES 
        source/cuda/meshtools.cu
        source/cuda/buoyancy.cu
        source/cuda/particle.cu
    )
endif()

# include dirs
SET(INCLUDE_PATHS
   source/pwrapper
   source/cuda
   source/util
   ${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source   
)

# Multithreading
if(MT)
    add_definitions( -DMT)    
    if(TBB)
        # Intel TBB
        add_definitions( -DTBB)

        if (DEBUG)
            add_definitions( -DTBB_USE_DEBUG=1 )
            list(APPEND F_LIBS tbb)
        else()
            list(APPEND F_LIBS tbb)
        endif()
        if (WIN32)
            find_package(TBB REQUIRED)
            list(APPEND INCLUDE_PATHS ${TBB_INCLUDE_DIRS})    
            list(APPEND F_LIB_PATHS ${TBB_LIBRARY_DIRS})
        elseif(APPLE)
            find_package(TBB REQUIRED)
            list(APPEND INCLUDE_PATHS ${TBB_INCLUDE_DIRS})    
            list(APPEND F_LIB_PATHS ${TBB_LIBRARY_DIRS})
        endif()
    else()
        # OpenMP
        add_definitions( -DOPENMP )
        
        if (WIN32)
            add_definitions( /openmp)
        else()
            add_definitions(-fopenmp)
            SET(F_LINKADD "${F_LINKADD} -fopenmp ")
        endif()
    endif()
endif()

#******************************************************************************
# add a target to generate API documentation with Doxygen

find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DX_PATH "doxy")
    foreach(it ${PP_SOURCES} ${PP_HEADERS} ${NOPP_SOURCES} ${NOPP_HEADERS})
        get_filename_component(CURPATH ${it} PATH)
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${DX_PATH}/${CURPATH}")
        set(CURDX "${DX_PATH}/${it}")
        string(REPLACE "/" "_" TGT ${CURDX})
        string(REPLACE "source/" "" INFILE ${it})
        add_custom_command(OUTPUT ${TGT} 
                        COMMAND prep docgen "0" ${MT_TYPE} "${CMAKE_CURRENT_SOURCE_DIR}/source/" "${INFILE}" "${CURDX}"
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

# Python
find_package(PythonLibs QUIET)
if((NOT PYTHONLIBS_FOUND) AND WIN32)
	set(PYTHON_INCLUDE_DIR "${WIN_PYTHON_PATH}/include")
	set(PYTHON_LIBRARY "${WIN_PYTHON_PATH}/libs/python${WIN_PYTHON_VERSION}.lib")	
endif()
find_package(PythonLibs REQUIRED)
list(APPEND INCLUDE_PATHS ${PYTHON_INCLUDE_DIRS})
list(APPEND F_LIBS ${PYTHON_LIBRARIES})

# Z compression
if (1)
	# default: build from source
	set(ZLIB_SRC adler32.c compress.c crc32.c deflate.c gzclose.c gzlib.c gzread.c gzwrite.c
				 inflate.c infback.c inftrees.c inffast.c trees.c uncompr.c zutil.c)
	foreach(it ${ZLIB_SRC})
        list(APPEND SILENT_SOURCES dependencies/zlib-1.2.8/${it})
    endforeach(it)
	list(APPEND INCLUDE_PATHS dependencies/zlib-1.2.8)
else()
	include(FindZLIB)
	list(APPEND INCLUDE_PATHS ${ZLIB_INCLUDE_DIR})
	list(APPEND F_LIBS ${ZLIB_LIBRARIES})
endif()

# CUDA
if(CUDA)
    add_definitions( -DCUDA=1 )
    find_package(CUDA QUIET REQUIRED)
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
    # ? if (USE64) set(CUDA_64_BIT_DEVICE_CODE ON) endif()
    if(WIN32)
        if(DEBUG)
            set(CUDA_NVCC_FLAGS -DDEBUG;--compiler-options;-Wall)
        else()
            set(CUDA_NVCC_FLAGS --use_fast_math;-DNDEBUG;--compiler-options;-Wall;-O2)
        endif()
    else()
        # CUDA does not support gcc > 4.5 yet
        if(DEBUG)
            set(CUDA_NVCC_FLAGS --pre-include;/usr/local/include/undef_atomics_int128.h;-g;-DDEBUG;-keep;--maxrregcount=31;--compiler-options;-Wall)
        else()
             #set(CUDA_NVCC_FLAGS -ccbin;gcc-4.5;--use_fast_math;-arch=sm_20;-DNDEBUG;--compiler-options;-Wall;-O3)       
             set(CUDA_NVCC_FLAGS --pre-include;/usr/local/include/undef_atomics_int128.h;--use_fast_math;-DNDEBUG;--compiler-options;-Wall;-O3)       
        endif()    
    endif()  
endif()

# increase FP precision?
if(DOUBLEPRECISION)
    add_definitions(-DFLOATINGPOINT_PRECISION=2)
endif()


#******************************************************************************
# generate repository info

# TEST! NT_DEBUG 
set(HGINFO "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source/hginfo.h")
MESSAGE(STATUS "Hg info target header ${HGINFO}")
add_custom_command(OUTPUT ${HGINFO}
				   COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tools/makeHgVersion.py" "${HGINFO}"
				   DEPENDS ${PP_SOURCES} ${PP_HEADERS} ${NOPP_SOURCES} ${NOPP_HEADERS} ${QT_SOURCES} ${QT_HEADERS}
				   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} )

#******************************************************************************
# apply preprocessor

set(SOURCES ${NOPP_SOURCES} ${SILENT_SOURCES})
set(HEADERS ${NOPP_HEADERS})
set(PP_REGCPP)
set(PP_REGS)
set(PP_PREPD "0")
if (PREPDEBUG)
    set(PP_PREPD "1")
endif()
foreach(it ${PP_SOURCES} ${PP_HEADERS})
    get_filename_component(CURPATH ${it} PATH)
    get_filename_component(CUREXT ${it} EXT)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${CURPATH}")
    set(CURPP "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${it}")
    string(REPLACE "source/" "" INFILE ${it})
    # preprocessor
    add_custom_command(OUTPUT ${CURPP}
                       COMMAND prep generate ${PP_PREPD} ${MT_TYPE} "${CMAKE_CURRENT_SOURCE_DIR}/source/" "${INFILE}" "${CURPP}"
                       DEPENDS prep
                       IMPLICIT_DEPENDS CXX ${it}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    if ("${CUREXT}" STREQUAL ".h" OR "${CUREXT}" STREQUAL ".py")
        list(APPEND PP_REGS "${CURPP}.reg")
        list(APPEND PP_REGCPP "${CURPP}.reg.cpp")
    endif()
    list(APPEND SOURCES ${CURPP})
endforeach(it)

# link reg files
add_custom_command(OUTPUT ${PP_REGCPP}
                  COMMAND prep link ${PP_REGS}
                  DEPENDS prep ${PP_REGS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Linking reg files")
#set_source_files_properties(${PP_REGS} ${PP_REGCPP} PROPERTIES GENERATED 1)
list(APPEND SOURCES ${PP_REGCPP})
    
#******************************************************************************
# QT for GUI

if(GUI)
    # remap
    set(QT_REMAP)
    foreach(it ${QT_HEADERS})
        list(APPEND QT_REMAP "${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/${it}")
    endforeach(it)
    add_definitions(-DGUI=1)
    list(APPEND INCLUDE_PATHS ${CMAKE_CURRENT_BINARY_DIR}
                ${CMAKE_CURRENT_BINARY_DIR}/${PP_PATH}/source/gui source/gui)
    
	cmake_policy(SET CMP0020 NEW)
    find_package(Qt5Core QUIET)
	if (Qt5Core_FOUND)
		message("Using Qt5")
        find_package(Qt5Widgets REQUIRED)
        find_package(Qt5OpenGL REQUIRED)
    
        qt5_wrap_cpp(MOC_OUTFILES ${QT_REMAP} )
        qt5_add_resources(QT_RES resources/res.qrc )
        add_definitions(-DGUI=1)
        add_definitions(${Qt5Widgets_DEFINITIONS})
        list(APPEND INCLUDE_PATHS ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS})
        list(APPEND F_LIBS ${Qt5Widgets_LIBRARIES} ${Qt5OpenGL_LIBRARIES})
        list(APPEND SOURCES ${MOC_OUTFILES} ${QT_RES})		
    else()
	    message("Using Qt4")
        find_package(Qt4 REQUIRED)
		set(QT_USE_QTOPENGL TRUE)

        qt4_wrap_cpp(MOC_OUTFILES ${QT_REMAP} )
        qt4_add_resources(QT_RES resources/res.qrc )
        include(${QT_USE_FILE})
        add_definitions(${QT_DEFINITIONS})
        list(APPEND F_LIBS ${QT_LIBRARIES})
        list(APPEND SOURCES ${MOC_OUTFILES} ${QT_RES})
    endif()

    if (APPLE)
		# mac opengl framework
        SET(F_LINKADD "${F_LINKADD} -framework OpenGL ")
    else()
		find_package(OpenGL REQUIRED)
		list(APPEND F_LIBS ${OPENGL_LIBRARIES})    
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
        add_definitions( -O3 -Wall -Wno-reorder)
    ENDIF()
ELSE()
    # get rid of some MSVC warnings
    add_definitions( /wd4018 /wd4146 /wd4800 )
	
	# for zlib
	add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
    add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
    
    # unsigned to signed conversions
    add_definitions( /wd4267 )
    # double <> single precision
    add_definitions( /wd4244 /wd4305 )
    # disable warnings for unsecure functions 
    add_definitions( /D "_CRT_SECURE_NO_WARNINGS" )	
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
	if (WIN32)
		# make nice folders for Visual Studio
		set_source_files_properties(${PP_SOURCES} ${PP_HEADERS} ${NOPP_HEADERS} PROPERTIES HEADER_FILE_ONLY TRUE) 
		add_executable(${EXECCMD} ${SOURCES} ${PP_REGISTRY} ${PP_SOURCES} ${PP_HEADERS})
		source_group(Generated FILES ${SOURCES} ${PP_REGISTRY} ${HEADERS})		
	else()
    	add_executable(${EXECCMD} ${SOURCES} ${PP_REGISTRY} ${HGINFO})
	endif()
	target_link_libraries( ${EXECCMD} ${F_LIBS} )
endif()

if(WIN32)
	# add no-build target to create nice MSVC folder
	add_custom_target(SourceFiles SOURCES ${PP_SOURCES} ${PP_HEADERS} ${QT_SOURCES} ${QT_HEADERS} ${NOPP_HEADERS} ${NOPP_SOURCES})
	source_group(Headers FILES ${PP_HEADERS} ${NOPP_HEADERS} ${QT_HEADERS})
	source_group(Sources FILES ${PP_SOURCES} ${NOPP_SOURCES} ${QT_SOURCES})	
endif()

# debug
MESSAGE(STATUS "DEBUG Flag-Summary - Includes: '${INCLUDE_PATHS}' | Libs: '${F_LIBS}' | LibPaths: '${F_LIB_PATHS}' ")


