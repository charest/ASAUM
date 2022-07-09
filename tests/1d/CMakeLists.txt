name_from_path(${CMAKE_CURRENT_SOURCE_DIR} _name)

# structured
set(_n 4)
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
  ARGS -i ${CMAKE_CURRENT_SOURCE_DIR}/sod.lua -p ${_n}
  PROCS 1 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
)

# structured with imex
set(_n 4)
set(_compare)
set(_standard)
append_file_list("csvfiles/000001/sod-imex.mesh.csv.%" ${_n} _compare)
append_file_list("csvfiles/000001/sod-imex.hydro.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000001/sod-imex.%.vts" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-imex.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-imex.hydro.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/vtkfiles/000001/sod-imex.%.vts" ${_n} _standard)

create_regression(
  NAME ${_name}-imex
  ARGS -i ${CMAKE_CURRENT_SOURCE_DIR}/sod-imex.lua -p ${_n}
  PROCS 1 2 3 4 5
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
  PROCS 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
)

# unstructured with naive partitioning
set(_compare)
set(_standard)
append_file_list("csvfiles/000001/sod-unstruct-naive.mesh.csv.%" ${_n} _compare)
append_file_list("csvfiles/000001/sod-unstruct-naive.hydro.csv.%" ${_n} _compare)
append_file_list("vtkfiles/000001/sod-unstruct-naive.%.vtu" ${_n} _compare)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct-naive.mesh.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/csvfiles/000001/sod-unstruct-naive.hydro.csv.%" ${_n} _standard)
append_file_list("${CMAKE_CURRENT_SOURCE_DIR}/vtkfiles/000001/sod-unstruct-naive.%.vtu" ${_n} _standard)

create_regression(
  NAME ${_name}-unstruct-naive
  ARGS 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/sod-unstruct.lua
    -p ${_n}
    -m ${CMAKE_CURRENT_SOURCE_DIR}/sod.g
    -a naive
    -c sod-unstruct-naive
  PROCS 1 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
)

# unstructured with pre-partitioned files
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
  PROCS 2 3 4 5
  COMPARE ${_compare}
  STANDARD ${_standard}
)

