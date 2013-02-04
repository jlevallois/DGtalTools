SET(DGTAL_TOOLS_SRC
  2dLocalEstimators
  3dLocalEstimators
  lengthEstimators
  vol2normalField
  curvatureBC
  curvatureMCMS
  tangentBC
  statisticsEstimators
  scaleSpaceIntegralInvariant2D
  scaleSpaceIntegralInvariant2D_ellipse
  scaleSpaceIntegralInvariant2D_flower
  scaleSpaceIntegralInvariant3D
)

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} )
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
