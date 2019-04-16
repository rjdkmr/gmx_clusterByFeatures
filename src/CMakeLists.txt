
add_library(gmx_clusterByFeatures MODULE pywrapper.cpp
                                         clusterbyfeatures.cpp
                                         logstream.cpp
                                         do_cluster.cpp
                                         distMat.cpp
                                         hole.cpp
                                         parseData.cpp
                                         )
target_link_libraries(gmx_clusterByFeatures PUBLIC ${GROMACS_LIBRARIES} PRIVATE pybind11::module)

if (INPLACE)
    set_target_properties(gmx_clusterByFeatures PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}"
                                        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/gmx_clusterByFeatures"
                                        LINKER_LANGUAGE CXX INSTALL_RPATH_USE_LINK_PATH ON)
else (INPLACE)
    set_target_properties(gmx_clusterByFeatures PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}"
                                        LINKER_LANGUAGE CXX INSTALL_RPATH_USE_LINK_PATH ON)
    INSTALL(TARGETS gmx_clusterByFeatures LIBRARY DESTINATION "${PYTHON_INSTDIR}/gmx_clusterByFeatures" COMPONENT python)
endif (INPLACE)
                                        
                                        


       
