SET(DGTAL_TOOLS_SRC
  2dLocalEstimators
  2dLocalEstimatorsII
  3dLocalEstimators
  2dLocalEstimators_0memory
  3dLocalEstimators_0memory
  3dLocalEstimators_0memory_old
  lengthEstimators
  vol2normalField
  curvatureBC
  curvatureMCMS
  tangentBC
  statisticsEstimators
  statisticsEstimators_0memory
)

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
