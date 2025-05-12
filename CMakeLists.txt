cmake_minimum_required(VERSION 3.10)
project(molmodel C CXX)

# Standards & warnings
set(CMAKE_C_STANDARD   99)
set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall)

# Dependencies
find_package(BLAS REQUIRED)       # Accelerate on macOS, OpenBLAS/BLAS on Linux
find_package(GSL  REQUIRED)       # GSL and its CBLAS

# FFTW3 via pkg-config if available, else manual probe
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_FFTW3 REQUIRED fftw3)
  add_library(FFTW3::fftw3 INTERFACE IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PC_FFTW3_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES        "${PC_FFTW3_LIBRARIES}"
  )
else()
  find_path(FFTW3_INCLUDE_DIR fftw3.h HINTS /usr/local/include /opt/homebrew/include)
  find_library(FFTW3_LIB fftw3 HINTS /usr/local/lib /opt/homebrew/lib)
  add_library(FFTW3::fftw3 INTERFACE IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES       "${FFTW3_LIB}"
  )
endif()

# Consolidated third-party libs
set(THIRD_PARTY
  ${BLAS_LIBRARIES}
  FFTW3::fftw3
  GSL::gsl
  GSL::gslcblas
  m
)

# Gather all non-driver sources
file(GLOB_RECURSE SRCS src/*.c src/*.cpp)
list(REMOVE_ITEM SRCS
  src/geoprog.cpp src/msimprog.cpp src/cgmin.c src/pppm-example.c
  src/rmsdprog.cpp src/histoprog.cpp
)

# Generate IO-based sources
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/gen)
set(IO_HDRS
  include/atom.hpp     include/msys.hpp    include/mdyn.hpp
  include/eleclj.hpp   include/checkgr.hpp include/callback.hpp
  include/rdf.hpp      include/intra.hpp   include/mmin.hpp
  include/msimprog.hpp include/geohist.hpp  include/acid.hpp
  include/constraint.hpp include/covalent.hpp include/contact.hpp
  include/rmsdcb.hpp   include/waterintra.hpp
  include/cmmsd.hpp    include/harmonic.hpp include/zmatrix.hpp
  include/espfit.hpp   include/bcifit.hpp   include/avgmf.hpp
  include/cmpole.hpp   include/ppot.hpp     include/hist.hpp
  include/prog.hpp
)
foreach(h IN LISTS IO_HDRS)
  get_filename_component(n ${h} NAME_WE)
  set(gen ${CMAKE_BINARY_DIR}/gen/${n}io.cpp)
  add_custom_command(
    OUTPUT            ${gen}
    COMMAND           ${CMAKE_SOURCE_DIR}/io -o ${h} > ${gen}
    DEPENDS           ${CMAKE_SOURCE_DIR}/io ${CMAKE_SOURCE_DIR}/${h}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT           "Generating ${n}io.cpp"
  )
  list(APPEND SRCS ${gen})
endforeach()

# Common library
add_library(molcommon STATIC ${SRCS})
target_include_directories(molcommon PUBLIC include)
target_link_libraries(molcommon PUBLIC ${THIRD_PARTY})

# Executables
add_executable(geo           src/geoprog.cpp)
add_executable(msim          src/msimprog.cpp)
add_executable(cgmin-example src/cgmin.c)
add_executable(pppm-example  src/pppm-example.c)
add_executable(rmsd          src/rmsdprog.cpp)
add_executable(hist          src/histoprog.cpp)
target_compile_definitions(cgmin-example PRIVATE SIMPLE_EXAMPLE)
target_compile_definitions(pppm-example PRIVATE SIMPLE_EXAMPLE)
foreach(e IN ITEMS geo msim cgmin-example pppm-example rmsd hist)
  target_link_libraries(${e} PRIVATE molcommon)
endforeach()

# Tests
enable_testing()
add_test(cgmin   cgmin-example)
add_test(pppm    pppm-example)
add_test(rmsd    rmsd ${CMAKE_SOURCE_DIR}/dat/v1 ${CMAKE_SOURCE_DIR}/dat/v2)
