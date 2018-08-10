add_executable(gmx_clusterByFeatures  gmx_clusterbyfeatures.h python_helper.cpp python_helper.h do_cluster.h do_cluster.cpp logstream.h logstream.cpp gmx_clusterbyfeatures.cpp main.cpp)
target_link_libraries(gmx_clusterByFeatures ${PYTHON_LIBRARIES} ${GROMACS_LIBRARIES} )
SET_TARGET_PROPERTIES(gmx_clusterByFeatures PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS gmx_clusterByFeatures DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
