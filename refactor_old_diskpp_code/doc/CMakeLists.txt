option(BUILD_DOC "Build documentation" OFF)

if (BUILD_DOC)
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
    message("Doxygen not found, I will not generate/install the documentation")
else()
    #message("Doxygen is found")

   configure_file(Doxyfile.in Doxyfile)

   set(DOXYGEN_INPUT ${CMAKE_BINARY_DIR}/doc/Doxyfile)
   set(DOXYGEN_OUTPUT ${APIDOC_DIR}/html/index.html)

   add_custom_target(diskpp_doc ALL
     COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT} > /dev/null
     COMMAND ${CMAKE_COMMAND} -E echo "Done."
     )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html
          DESTINATION share/doc/${CMAKE_PROJECT_NAME})

endif()
endif(BUILD_DOC)
