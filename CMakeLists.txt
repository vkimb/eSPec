cmake_minimum_required (VERSION 3.7 FATAL_ERROR)

project (fsp_v2)

enable_language (Fortran C)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE  "Release")
endif (NOT CMAKE_BUILD_TYPE)

#get_filename_component (default_prefix ".." ABSOLUTE)
  
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fpp -O3  -ipo -mkl=sequential -diag-disable=11021,8291,8293")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fpp -O0  -g -mkl=sequential -check all -diag-disable=11021,8291,8293")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message (" * Picking up a Debug build * ")
  #set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fpp -O0 -openmp -g -mkl -check all -diag-disable=11021,8291,8293 ")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  message (" * Picking up a Release build * ")
  #set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}   -fpp -O3 -openmp -ipo -mkl -diag-disable=11021,8291,8293 ")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")
    
include_directories ( /opt/sundials/include  /opt/intel/mkl/include  )
link_directories (  /opt/intel/mkl/lib/intel64 /opt/intel/mkl/lib  /opt/sundials/lib)
#include_directories ( /opt/hdf5/include /opt/openmpi_intel/include  /opt/intel/mkl/include  /opt/petsc_3.7.3/include  /opt/slepc_3.7.2/include )
#link_directories ( /opt/hdf5/lib /opt/openmpi_intel/lib /opt/intel/mkl/lib/intel64 /opt/intel/mkl/lib /opt/petsc_3.7.3/lib /opt/slepc_3.7.2/lib )

#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp -mkl")
#add_definitions(-mkl -fpp )
# use, i.e. don't skip the full RPATH for the build tree

add_executable(espec.x src/espec.f src/rdinput.f src/compar.f src/chlength.f src/rdpt.f src/gvpot.f src/potfuncs.f src/mtrxdiag.f src/lanczs.f src/lnz.f src/reort.f src/ecnorm.f src/au.f src/au1d.f src/au2di.f src/au2dit.f src/au2dil.f src/au2dc.f src/au3d.f src/rans.f src/lanczsc.f src/rdpte.f src/dfdxi.f src/rif.f src/sod.f src/psod.f src/sil.f src/plnz.f src/fcorrel.f src/eigenerg.f src/spectrumtd.f src/spectrumti.f src/transdipole.f src/makecst.f src/getcorr.f src/prtpot.f src/prteigvc.f src/ef.f src/chpot.f src/ppsod.f src/pplnz.f src/hwhm.f src/dipl.f src/prpt2.f  src/rdpt2.f src/sodfft.f src/psodfft.f src/ddmtv.f src/ddmtvtn.f src/dvtn.f src/spofft.f src/pspofft.f src/dkeft.f src/eigenergfft.f src/eigenergfft1.f src/absorb.f src/s2ppsod.f src/abm1.f src/s2ppabm.f src/sod1.f src/s2ppabm2.f src/init_cond.f src/gauss.f src/ap1d.f src/ap2d.f src/ap2dct.f src/au2dct.f src/coupled3.f90 src/rdpt-min.f src/rhside.f90 src/csection.c src/splinesurf.f src/io.f90 src/rtlib.f90 src/cfft.f src/fft3d.f )
target_link_libraries(espec.x iomp5  fftw3 mkl_sequential mkl_blas95_ilp64 mkl_lapack95_ilp64 mkl_intel_ilp64 sundials_fcvode sundials_fnvecserial sundials_cvode sundials_nvecserial)
