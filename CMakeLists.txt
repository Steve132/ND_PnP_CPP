cmake_minimum_required(VERSION 3.10)

project(PnPDynamic)

add_library(PnPDynamic
	PnPDynamic.hpp
	PnPDynamic.cpp
)
target_include_directories(PnPDynamic PUBLIC deps)
set_target_properties(PnPDynamic PROPERTIES CXX_STANDARD 11)
if(NOT MSVC)
#target_compile_options(PnPDynamic PUBLIC "-march=native")
endif()
add_executable(pnp_tester
	pnp_tester.cpp
)

set_target_properties(pnp_tester PROPERTIES CXX_STANDARD 11)
target_include_directories(pnp_tester PRIVATE deps)
target_link_libraries(pnp_tester PnPDynamic)
