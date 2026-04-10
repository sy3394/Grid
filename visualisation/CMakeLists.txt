cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(GridViewer)

#list(APPEND CMAKE_PREFIX_PATH "/Users/peterboyle/QCD/vtk/VTK-9.4.2-install/")
list(APPEND CMAKE_PREFIX_PATH "/ccs/home/syamamoto/src/vtk/build/")

find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonTransforms
  FiltersCore
  FiltersModeling
  IOImage
  IOFFMPEG
  InteractionStyle
  InteractionWidgets
  RenderingAnnotation
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "GridViewer: Unable to find the VTK build folder.")
endif()

# HDF5 (C library only — C++ bindings may be absent on some systems).
# CMake's find_package may not locate it on Cray PE; fall back to env vars
# set by "module load cray-hdf5" (CRAY_HDF5_DIR, HDF5_DIR, HDF5_ROOT).
find_package(HDF5 COMPONENTS C QUIET)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5 found via CMake: ${HDF5_INCLUDE_DIRS}")
else()
  foreach(_hdf5_env CRAY_HDF5_DIR HDF5_DIR HDF5_ROOT)
    if (DEFINED ENV{${_hdf5_env}})
      set(_hdf5_inc "$ENV{${_hdf5_env}}/include")
      if (EXISTS "${_hdf5_inc}/hdf5.h")
        include_directories(${_hdf5_inc})
        message(STATUS "HDF5 include from env ${_hdf5_env}: ${_hdf5_inc}")
        break()
      endif()
    endif()
  endforeach()
endif()


# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(FieldDensityAnimate MACOSX_BUNDLE FieldDensityAnimate.cxx )
target_link_libraries(FieldDensityAnimate PRIVATE ${VTK_LIBRARIES})

add_executable(FieldDensityAnimateMultiFiles MACOSX_BUNDLE FieldDensityAnimateMultiFiles.cxx )
target_link_libraries(FieldDensityAnimateMultiFiles PRIVATE ${VTK_LIBRARIES})

add_executable(FieldDensityFindRegion MACOSX_BUNDLE FieldDensityFindRegion.cxx )
target_link_libraries(FieldDensityFindRegion PRIVATE ${VTK_LIBRARIES})

add_executable(FieldDensityEigen MACOSX_BUNDLE FieldDensityEigen.cxx )
target_link_libraries(FieldDensityEigen PRIVATE ${VTK_LIBRARIES})

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS FieldDensityAnimate FieldDensityAnimateMultiFiles FieldDensityFindRegion
  MODULES ${VTK_LIBRARIES}
)
