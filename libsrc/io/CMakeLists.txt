target_sources( asaum PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/csv_writer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/file_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/side_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_writer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/write_binary.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/write_zlib.cpp )

if (HAVE_EXODUS)
  target_sources( asaum PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/exo_reader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exo_types.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exo_writer.cpp
  )
endif()
