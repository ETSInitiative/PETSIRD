find_package(xtensor REQUIRED)
# find_package(xtensor-blas REQUIRED)
find_package(BLAS REQUIRED)

# create a petsird_helpers library target
# Currently it is just empty, but allows creating target properties which are transitive
add_library(petsird_helpers INTERFACE)

target_include_directories(petsird_helpers
  INTERFACE
  ${xtensor_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(petsird_helpers INTERFACE petsird BLAS::BLAS)

add_executable(petsird_generator petsird_generator.cpp)
target_link_libraries(petsird_generator PRIVATE petsird_helpers)

add_executable(petsird_analysis petsird_analysis.cpp)
target_link_libraries(petsird_analysis PRIVATE petsird_helpers)

foreach(t petsird petsird_generated petsird_helpers)
  if(TARGET ${t})
    get_target_property(_libs ${t} INTERFACE_LINK_LIBRARIES)
    if(_libs)
      list(FILTER _libs EXCLUDE REGEX "xtensor")
      set_target_properties(${t}
        PROPERTIES INTERFACE_LINK_LIBRARIES "${_libs}")
    endif()
  endif()
endforeach()
