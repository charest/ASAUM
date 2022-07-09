
add_executable( example main.cpp driver.cpp )
target_link_libraries(example asaum)

set_property(TARGET example PROPERTY CXX_STANDARD 14)
set_property(TARGET example PROPERTY CXX_STANDARD_REQUIRED on)
set_property(TARGET example PROPERTY CXX_EXTENSIONS off)

