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
# Try CMake's finder first; if it fails fall back to env vars and then to
# the known Frontier spack path baked into Grid's configure summary.
find_package(HDF5 COMPONENTS C QUIET)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5 found via CMake: ${HDF5_INCLUDE_DIRS}")
else()
  set(_hdf5_inc "")
  # Env vars set by "module load cray-hdf5" or spack
  foreach(_hdf5_env CRAY_HDF5_DIR HDF5_DIR HDF5_ROOT)
    if (DEFINED ENV{${_hdf5_env}} AND EXISTS "$ENV{${_hdf5_env}}/include/hdf5.h")
      set(_hdf5_inc "$ENV{${_hdf5_env}}/include")
      break()
    endif()
  endforeach()
  # Frontier spack fallback (from Grid's configure summary)
  if (NOT _hdf5_inc)
    set(_frontier_hdf5 "/autofs/nccs-svm1_home2/syamamoto/src/spack/opt/spack/linux-zen3/hdf5-1.14.6-zjapmrieixabrttqp6zmixtwzuxymign/include")
    if (EXISTS "${_frontier_hdf5}/hdf5.h")
      set(_hdf5_inc "${_frontier_hdf5}")
    endif()
  endif()
  if (_hdf5_inc)
    include_directories(${_hdf5_inc})
    # Derive lib dir from include dir (../lib)
    get_filename_component(_hdf5_root "${_hdf5_inc}" DIRECTORY)
    find_library(_hdf5_lib hdf5 HINTS "${_hdf5_root}/lib" NO_DEFAULT_PATH)
    if (_hdf5_lib)
      set(HDF5_LIBRARIES ${_hdf5_lib})
    else()
      set(HDF5_LIBRARIES "-L${_hdf5_root}/lib -lhdf5")
    endif()
    message(STATUS "HDF5 include: ${_hdf5_inc}")
    message(STATUS "HDF5 library: ${HDF5_LIBRARIES}")
  else()
    message(WARNING "HDF5 not found — HDF5 read path will be disabled")
  endif()
endif()


# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(FieldDensityAnimate MACOSX_BUNDLE FieldDensityAnimate.cxx )
target_link_libraries(FieldDensityAnimate PRIVATE ${VTK_LIBRARIES})

add_executable(FieldDensityAnimateMultiFiles MACOSX_BUNDLE FieldDensityAnimateMultiFiles.cxx )
target_link_libraries(FieldDensityAnimateMultiFiles PRIVATE ${VTK_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(FieldDensityFindRegion MACOSX_BUNDLE FieldDensityFindRegion.cxx )
target_link_libraries(FieldDensityFindRegion PRIVATE ${VTK_LIBRARIES})

add_executable(FieldDensityEigen MACOSX_BUNDLE FieldDensityEigen.cxx )
target_link_libraries(FieldDensityEigen PRIVATE ${VTK_LIBRARIES} ${HDF5_LIBRARIES})

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS FieldDensityAnimate FieldDensityAnimateMultiFiles FieldDensityFindRegion
  MODULES ${VTK_LIBRARIES}
)
