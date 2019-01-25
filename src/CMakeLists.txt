
add_library(gmx_clusterByFeatures MODULE gmx_clusterbyfeatures.h do_cluster.h do_cluster.cpp logstream.h logstream.cpp gmx_clusterbyfeatures.cpp pywrapper.cpp)
target_link_libraries(gmx_clusterByFeatures PUBLIC ${GROMACS_LIBRARIES} PRIVATE pybind11::module)
set_target_properties(gmx_clusterByFeatures PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}"
                                        LIBRARY_OUTPUT_DIRECTORY ${gmx_clusterByFeatures_PYTHON_OUTPUT_DIRECTORY} 
                                        LINKER_LANGUAGE CXX)
                                          
                                          
INSTALL(TARGETS gmx_clusterByFeatures
        LIBRARY DESTINATION ${gmx_clusterByFeatures_PythonDir}
        COMPONENT python
       )
