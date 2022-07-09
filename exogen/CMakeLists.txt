
add_executable( exogen exogen.cpp )
target_link_libraries(exogen asaum)

set_property(TARGET exogen PROPERTY CXX_STANDARD 14)
set_property(TARGET exogen PROPERTY CXX_STANDARD_REQUIRED on)
set_property(TARGET exogen PROPERTY CXX_EXTENSIONS off)

