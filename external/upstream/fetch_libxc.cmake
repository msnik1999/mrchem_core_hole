cpm_set_find_behaviour(${LIBXC_FIND_BEHAVIOUR})
set(_saved_install_includedir "${CMAKE_INSTALL_INCLUDEDIR}")
set(CMAKE_INSTALL_INCLUDEDIR "include/LibXC" CACHE STRING "" FORCE)
CPMAddPackage(
  NAME Libxc
  VERSION 7
  GIT_TAG 7.0.0
  GITLAB_REPOSITORY libxc/libxc
  FIND_PACKAGE_ARGUMENTS "CONFIG NO_CMAKE_PATH NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PACKAGE_REGISTRY"
  OPTIONS
  "BUILD_TESTING OFF"
  "ENABLE_TESTS OFF"
  "BUILD_SHARED_LIBS ON"
  "ENABLE_FORTRAN OFF"
  "DISABLE_FXC ON"
  "DISABLE_KXC ON"
  "DISABLE_LXC ON"
  )
set(CMAKE_INSTALL_INCLUDEDIR "${_saved_install_includedir}" CACHE STRING "" FORCE)

if(TARGET xc)
  if(NOT TARGET Libxc::xc)
    add_library(Libxc::xc ALIAS xc)
  endif()
  target_include_directories(xc
    INTERFACE
      $<BUILD_INTERFACE:${libxc_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}>
      $<INSTALL_INTERFACE:include/LibXC>
    )
endif()
