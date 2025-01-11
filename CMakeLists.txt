# much of this is owed to Leighton Wilson

cmake_minimum_required(VERSION 3.9...3.14)
project(FARRSIGHT LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the C++ compiler (e.g., NVIDIA's pgc++)
set(CMAKE_CXX_COMPILER nvc++)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release" "None")
endif()

set(default_include /usr/local/include)
if(NOT INCLUDE_HEADERS)
    message(STATUS "Setting include directory to '${default_include}' as none was specified.")
    set(INCLUDE_HEADERS ${default_include} CACHE FILEPATH "Choose include directory")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(ENABLE_GPU_BUILD "build OpenACC GPU version"             ON)

find_package(OpenMP     REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
# Interface for OpenACC
add_library(OpenACC_CXX INTERFACE)
set_property(TARGET OpenACC_CXX PROPERTY INTERFACE_COMPILE_OPTIONS -acc -Minfo=accel)


# from Leighton
#Some PGI versions seem to require a link flag for OpenMP
set_property(TARGET OpenMP::OpenMP_C APPEND
                                     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         $<$<CXX_COMPILER_ID:NVHPC>:-mp>) 
# #Creating an interface for OpenACC
# add_library(OpenACC_CXX INTERFACE)
# set_property(TARGET OpenACC_CXX     PROPERTY
#                                          INTERFACE_COMPILE_OPTIONS
#                                          -acc
# #                                         -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
#                                          -Minfo)
# set_property(TARGET OpenACC_CXX     PROPERTY
#                                          INTERFACE_LINK_LIBRARIES
#                                          -acc
# #                                         -ta=tesla:ccall,cc60
#                                         )
# if(ENABLE_GPU_BUILD AND NOT (CMAKE_C_COMPILER_ID MATCHES NVHPC))
#     message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
#     set(ENABLE_GPU_BUILD OFF CACHE BOOL "build OpenACC GPU version" FORCE)
# endif()


add_subdirectory(src)


# add the executable
add_executable(farrsight_cpu run_FARRSIGHT2.cpp)
target_link_libraries     (farrsight_cpu PUBLIC   OpenMP::OpenMP_CXX
                                                SimObj
                                                FieldObj
                                                PanelObj
                                                icsObj
                                                AMRObj_cpu)
target_include_directories(farrsight_cpu PUBLIC ${INCLUDE_HEADERS})
target_include_directories(farrsight_cpu PUBLIC   inc)


install(TARGETS farrsight_cpu  RUNTIME DESTINATION bin
                                LIBRARY DESTINATION lib
                            ARCHIVE DESTINATION lib)

# Find the Boost package
find_package(Boost REQUIRED)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

#TARGET_LINK_LIBRARIES(farrsight_cpu LINK_PUBLIC ${BOOST_LIB} )
TARGET_LINK_LIBRARIES(farrsight_cpu LINK_PUBLIC ${Boost_LIBRARIES} )




# Set additional flags for NVIDIA compilers if ENABLE_GPU_BUILD is on
if(ENABLE_GPU_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -acc -Minfo=accel -gpu=managed") # -gpu=managed -ta=tesla:managed -gpu=nordc
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -acc")
else()
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "Build OpenACC GPU version" FORCE)
endif()


if (ENABLE_GPU_BUILD)
    add_executable(farrsight_gpu run_FARRSIGHT2.cpp) 
    #find_library(BARYTREE_LIB_GPU BaryTree_gpu REQUIRED)
    target_link_libraries(farrsight_gpu PUBLIC   OpenMP::OpenMP_CXX
                                                OpenACC_CXX
                                                SimObj_gpu
                                                FieldObj_gpu
                                                PanelObj
                                                icsObj
                                                AMRObj_cpu)
    target_include_directories(farrsight_gpu PUBLIC ${INCLUDE_HEADERS})
    target_include_directories(farrsight_gpu PUBLIC   inc)
    target_compile_definitions(farrsight_gpu PUBLIC OPENACC_ENABLED)

    install(TARGETS farrsight_gpu RUNTIME DESTINATION bin
                                LIBRARY DESTINATION lib
                                ARCHIVE DESTINATION lib)

TARGET_LINK_LIBRARIES(farrsight_gpu LINK_PUBLIC ${Boost_LIBRARIES} )
                               
endif()



#enable_testing()

#add_test(NAME Runs COMMAND ${EXE} 1 0.0 1.0 -2.0 2.0)
