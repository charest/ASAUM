target_sources( asaum PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/bamr/bamr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bamr/bamr_block.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/box_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/file_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/block_ijk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/block_ijk_neighbor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/ctm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/mesh_ijk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/range.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/sector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ijk/transform.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_block.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_boundary.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_struct_block.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_struct_block_cpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_struct_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_unstruct_block.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_unstruct_block_cpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/partition.cpp
)
  
if (HAVE_EXODUS)
  target_sources( asaum PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/exodus_mesh.cpp
  )
endif()
