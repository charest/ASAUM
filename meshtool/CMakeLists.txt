
add_executable( meshtool main.cpp driver.cpp )
target_link_libraries(meshtool asaum)

set_property(TARGET meshtool PROPERTY CXX_STANDARD 14)
set_property(TARGET meshtool PROPERTY CXX_STANDARD_REQUIRED on)
set_property(TARGET meshtool PROPERTY CXX_EXTENSIONS off)

