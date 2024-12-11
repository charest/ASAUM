name_from_path(${CMAKE_CURRENT_SOURCE_DIR} _name)

# structured
set(_n 4)
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vts" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/struct/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/struct/vtkfiles/000000/out.%.vts" ${_n} _standard)

create_regression(
  NAME ${_name}-struct
  ARGS -i ${CMAKE_CURRENT_SOURCE_DIR}/struct.lua -p ${_n}
  PROCS 1 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY struct
)

# unstructured
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct/vtkfiles/000000/out.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct
  ARGS 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g
  PROCS 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct
)

# unstructured with naive partitioning
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-naive/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-naive/vtkfiles/000000/out.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct-naive
  ARGS 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g
    -a naive
  PROCS 1 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct-naive
)

# unstructured with pre-partitioned files
set(_compare)
set(_standard)
append_file_list("csvfiles/000000/out.mesh.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000000/out.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-prepart/csvfiles/000000/out.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/unstruct-prepart/vtkfiles/000000/out.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct-prepart
  ARGS 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/mesh.g.02
  PROCS 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
  WORKING_DIRECTORY unstruct-prepart
)

