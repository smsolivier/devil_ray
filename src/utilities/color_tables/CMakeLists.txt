set(color_tables_thirdparty_libs dray)

################################################
# openmp support
################################################
if(ENABLE_OPENMP)
  list(APPEND color_tables_thirdparty_libs openmp)
endif()

blt_add_executable(
  NAME color_tables
  SOURCES color_tables.cpp
  DEPENDS_ON ${color_tables_thirdparty_libs}
  OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
