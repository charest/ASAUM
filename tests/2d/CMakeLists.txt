name_from_path(${CMAKE_CURRENT_SOURCE_DIR} _name)

# structured
set(_n 9)
set(_compare)
set(_standard)
append_file_list("csvfiles/000001/sod.mesh.csv.%" ${_n} _compare)
append_file_list("csvfiles/000001/sod.hydro.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000001/sod.%.vts" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod.hydro.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/vtkfiles/000001/sod.%.vts" ${_n} _standard)

create_regression(
  NAME ${_name}
  ARGS -i ${CMAKE_CURRENT_SOURCE_DIR}/sod.lua -p 9
  PROCS 1 3 9 10
  COMPARE ${_compare}
  STANDARD ${_standard}
)

# unstructured
set(_compare)
set(_standard)
append_file_list("csvfiles/000001/sod-unstruct.mesh.csv.%" ${_n} _compare)
append_file_list("csvfiles/000001/sod-unstruct.hydro.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000001/sod-unstruct.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct.hydro.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/vtkfiles/000001/sod-unstruct.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct
  ARGS
    -i ${CMAKE_CURRENT_SOURCE_DIR}/sod-unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/sod.g
  PROCS 10
  COMPARE ${_compare}
  STANDARD ${_standard}
)

# unstructured - prepart
set(_compare)
set(_standard)
append_file_list("csvfiles/000001/sod-unstruct-prepart.mesh.csv.%" ${_n} _compare)
append_file_list("csvfiles/000001/sod-unstruct-prepart.hydro.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000001/sod-unstruct-prepart.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct-prepart.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct-prepart.hydro.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/vtkfiles/000001/sod-unstruct-prepart.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct-prepart
  ARGS
    -i ${CMAKE_CURRENT_SOURCE_DIR}/sod-unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/sod.g.02
    -c sod-unstruct-prepart
  PROCS 3 9 10
  COMPARE ${_compare}
  STANDARD ${_standard}
)

