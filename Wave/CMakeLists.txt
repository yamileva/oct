cmake_minimum_required(VERSION 3.10)

# Set the project name
project(CMakeWave)
set(CMAKE_CXX_STANDARD 20)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    #src/Cpython/Call_Cpython.cpp
    src/Solver/Solver.cpp
    #src/Utils/Converter.cpp
    src/Utils/Grid.cpp
    src/CMakeWave.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME} ${SOURCES})

option(ENABLE_PYTHON "Enable Python" OFF)
if(ENABLE_PYTHON)
    find_package(Python3 3.7 COMPONENTS Interpreter NumPy REQUIRED)
    include_directories(${PYTHON_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES} Python3::NumPy)
endif(ENABLE_PYTHON)

option(ENABLE_HDF5 "Enable HDF5" ON)
if(ENABLE_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    include_directories(${HDF5_INCLUDE_DIRS})
    set(CMAKE_C_FLAGS "${CMAKE_FLAGS} ${MPI_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
    target_link_libraries( ${PROJECT_NAME} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
endif(ENABLE_HDF5)


# Set the directories that should be included in the build command for this target
target_include_directories(${PROJECT_NAME}
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

#set(MSVC_COMPILE_RELEASE_OPTIONS "/Od")
#target_compile_options(CMakeWave PRIVATE  "$<$<CONFIG:RELEASE>:${MSVC_COMPILE_RELEASE_OPTIONS}>")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    set(CMAKE_VERBOSE_MAKEFILE 1)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_DEFINITIONS "Bits64_;UNIX;_BOOL;LINUX;FUNCPROTO;_GNU_SOURCE;LINUX_64;REQUIRE_IOSTREAM")

    set(GCC_COMPILE_OPTIONS "-m64;-fPIC;-fno-strict-aliasing;-Wall;-Wno-multichar;-Wno-comment;-Wno-sign-compare;-funsigned-char;-pthread;-Wno-deprecated;-Wno-reorder;-ftemplate-depth-64;-fno-gnu-keywords;-std=c++0x;-Winline")
    set(GCC_COMPILE_DEBUG_OPTIONS "${GCC_COMPILE_OPTIONS};-ggdb;-O0")
    set(GCC_COMPILE_RELEASE_OPTIONS "${GCC_COMPILE_OPTIONS};-O3")
    target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${GCC_COMPILE_DEBUG_OPTIONS}>")
    target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${GCC_COMPILE_RELEASE_OPTIONS}>")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
    set(MSVC_COMPILE_OPTIONS "/MP;/W3;/w34710;/Gy;/Zc:wchar_t;/nologo;/std:c++0x /EHsc")
    set(MSVC_COMPILE_DEBUG_OPTIONS "${MSVC_COMPILE_OPTIONS} /ZI /Od;/openmp")
    set(MSVC_COMPILE_RELEASE_OPTIONS "${MSVC_COMPILE_OPTIONS} /O2;/openmp")
    target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${MSVC_COMPILE_DEBUG_OPTIONS}>")
    target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${MSVC_COMPILE_RELEASE_OPTIONS}>")
endif()


#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


#message(STATUS "CMAKE_CXX_FLAGS_RELEASE =${CMAKE_CXX_FLAGS_RELEASE}" )