include_directories(${PROJECT_SOURCE_DIR}/inc)

add_library(icsObj initial_distributions.cpp)
target_link_libraries(icsObj      PUBLIC  $<$<CXX_COMPILER_ID:NVIDIA>:m>) 
target_include_directories(icsObj PUBLIC ${INCLUDE_HEADERS} ${Boost_INCLUDE})

# Find the Boost package
find_package(Boost REQUIRED)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

#TARGET_LINK_LIBRARIES(farrsight_cpu LINK_PUBLIC ${BOOST_LIB} )
TARGET_LINK_LIBRARIES(icsObj LINK_PUBLIC ${Boost_LIBRARIES} )



add_library(PanelObj Panel.cpp)

add_library(AMRObj_cpu  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    AMRStructure_interp.cpp)
target_link_libraries(AMRObj_cpu    PUBLIC  PanelObj
                                            icsObj
                                            OpenMP::OpenMP_CXX
                                            $<$<CXX_COMPILER_ID:NVIDIA>:m>)
target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
                                                ${INCLUDE_HEADERS})