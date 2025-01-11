include_directories(${PROJECT_SOURCE_DIR}/inc)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Interface for OpenACC
# add_library(OpenACC_CXX INTERFACE)
# set_property(TARGET OpenACC_CXX PROPERTY INTERFACE_COMPILE_OPTIONS -acc -Minfo=accel)


add_library(FieldObj FieldStructure.cpp)
target_link_libraries(FieldObj       PUBLIC 
                                            OpenMP::OpenMP_CXX
                                            $<$<CXX_COMPILER_ID:NVIDIA>:m>)
target_include_directories(FieldObj PUBLIC ${INCLUDE_HEADERS})

# Set additional flags for NVIDIA compilers if ENABLE_GPU_BUILD is on
if(ENABLE_GPU_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -acc -Minfo=accel -gpu=managed") # -gpu=managed -ta=tesla:managed -gpu=nordc
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -acc")
else()
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "Build OpenACC GPU version" FORCE)
endif()

if(ENABLE_GPU_BUILD)
    add_library(FieldObj_gpu FieldStructure.cpp)
    target_compile_definitions(FieldObj_gpu PUBLIC OPENACC_ENABLED)
    target_link_libraries(FieldObj_gpu  PUBLIC  
                                                OpenMP::OpenMP_CXX
                                                OpenACC_CXX
                                                $<$<CXX_COMPILER_ID:NVIDIA>:m>)
    # target_compile_options(FieldObj_gpu PRIVATE -O0)
    target_include_directories(FieldObj_gpu PUBLIC ${INCLUDE_HEADERS} inc)
endif()

add_subdirectory(structure_files)
add_subdirectory(simulation_files)

