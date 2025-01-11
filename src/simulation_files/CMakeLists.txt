include_directories(${PROJECT_SOURCE_DIR}/inc)



add_library(SimObj AMRSimulation.cpp
                   AMRSimulation_load.cpp
                   AMRSimulation_gather_scatter.cpp
                   AMRSimulation_field.cpp
                   AMRSimulation_step.cpp
                   AMRSimulation_remesh.cpp
                   AMRSimulation_run.cpp
                   AMRSimulation_io.cpp)
target_link_libraries(SimObj       PUBLIC  AMRObj_cpu
                                            FieldObj
                                            OpenMP::OpenMP_CXX
                                            $<$<CXX_COMPILER_ID:NVIDIA>:m>)
target_include_directories(SimObj PUBLIC ${INCLUDE_HEADERS} ${Boost_INCLUDE})

# Find the Boost package
find_package(Boost REQUIRED)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(SimObj LINK_PUBLIC ${Boost_LIBRARIES} )


if (ENABLE_GPU_BUILD)

    add_library(SimObj_gpu AMRSimulation.cpp
                    AMRSimulation_load.cpp
                    AMRSimulation_gather_scatter.cpp
                    AMRSimulation_field.cpp
                    AMRSimulation_step.cpp
                    AMRSimulation_remesh.cpp
                    AMRSimulation_run.cpp
                    AMRSimulation_io.cpp)
    target_compile_definitions(SimObj_gpu PRIVATE OPENACC_ENABLED)
    target_link_libraries(SimObj_gpu       PUBLIC  AMRObj_cpu
                                                FieldObj_gpu
                                                OpenMP::OpenMP_CXX
                                                OpenACC_CXX
                                                $<$<CXX_COMPILER_ID:NVIDIA>:m>)
    target_include_directories(SimObj_gpu PUBLIC ${INCLUDE_HEADERS})

endif()