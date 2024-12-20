name_from_path(${CMAKE_CURRENT_SOURCE_DIR} _name)

# structured
set(_n 9)
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vts" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/struct/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/struct/vtkfiles/000000/out.%.vts" ${_n} _standard)

create_regression(
  NAME ${_name}-struct
  ARGS -i ${CMAKE_CURRENT_SOURCE_DIR}/struct.lua -p 9
  PROCS 1 3 9 10
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY struct
)

# unstructured
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/p1/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/p1/vtkfiles/000000/out.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct
  ARGS
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g
  PROCS 1
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct
)
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/p10/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/p10/vtkfiles/000000/out.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct
  ARGS
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g
  PROCS 10
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct
)

# unstructured - prepart
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-prepart/csvfiles/000000/sod-unstruct-prepart.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-prepart/vtkfiles/000000/sod-unstruct-prepart.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct-prepart
  ARGS
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g.02
  PROCS 3 9 10
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct-prepart
)

