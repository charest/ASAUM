set(LUA_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/fileutil.lua
  ${CMAKE_CURRENT_SOURCE_DIR}/interp.lua
  ${CMAKE_CURRENT_SOURCE_DIR}/stringutil.lua
)

foreach(_file ${LUA_FILES})
  get_filename_component(_name ${_file} NAME_WLE)
  set(LUA_DEF "${LUA_DEF} {\"${_name}\", \"${_file}\"},\n")
endforeach()

configure_file(lua_modules.hpp.in lua_modules.hpp @ONLY)
target_include_directories( asaum PRIVATE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> )

target_sources( asaum PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/lua.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lua_base.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lua_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lua_utils.cpp
)
