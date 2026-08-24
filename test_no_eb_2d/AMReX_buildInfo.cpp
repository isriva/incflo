
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2026-08-24 16:29:06.788949";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/data/home/mhudes/Code/incflo_ishan/test_no_eb_2d";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux gigan 6.8.0-84-generic #84-Ubuntu SMP PREEMPT_DYNAMIC Fri Sep  5 22:36:38 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "../../amrex_john";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "13.3.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "mpicxx";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "mpif90";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "-Werror=return-type -g1 -O3 -finline-limit=43210 -std=c++20  -pthread   -DBL_USE_MPI -DAMREX_USE_MPI -DBL_NO_FORT -DAMREX_GPU_MAX_THREADS=0 -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DAMREX_USE_FFT -DAMREX_XSDK -DNDEBUG -DOMPI_SKIP_MPICXX -DINCFLO_USE_FFT -Itmp_build_dir/s/2d.gnu.MPI.EXE -I. -I../../amrex_john/Src/FFT -I../src -I../src/boundary_conditions -I../src/convection -I../src/derive -I../src/diffusion -I../src/prob -I../src/projection -I../src/rheology -I../src/setup -I../src/utilities -I../src/analysis -I../../AMReX-Hydro/BDS -I../../AMReX-Hydro/Godunov -I../../AMReX-Hydro/MOL -I../../AMReX-Hydro/Projections -I../../AMReX-Hydro/Utils -I../../amrex_john/Src/Base -I../../amrex_john/Src/Base/Parser -I../../amrex_john/Src/AmrCore -I../../amrex_john/Src/Boundary -I../../amrex_john/Src/LinearSolvers/MLMG -I../../amrex_john/Tools/C_scripts ";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "-g1 -O3 -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none ";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "-L. ";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = " -lfftw3f -lfftw3 -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_cxx -lmpi ";
  return libraries;
}

const char* buildInfoGetMakeFlags() {

  static const char make_flags[] = "";
  return make_flags;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int const num_modules = X;
  int const num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "24.11-121-gedc5746-dirty";
  static const char HASH2[] = "20.04-3765-g2e25989a4";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

#ifdef AMREX_USE_CUDA
const char* buildInfoGetCUDAVersion() {

  static const char CUDA_VERSION[] = "";
  return CUDA_VERSION;
}
#endif

}
