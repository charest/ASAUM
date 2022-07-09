target_include_directories( asaum PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> )

#------------------------------------------------------------------------------#
# configure header
#------------------------------------------------------------------------------#

configure_file(config.hpp.in config.hpp @ONLY)

#------------------------------------------------------------------------------#
# add source subdirectories
#------------------------------------------------------------------------------#
add_subdirectory( comm )
add_subdirectory( common )
add_subdirectory( io )
add_subdirectory( lua )
add_subdirectory( math )
add_subdirectory( mesh )
add_subdirectory( utils )

