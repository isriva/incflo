#include <incflo.H>
#include <AMReX_buildInfo.H>
#include "chrono"

void writeBuildInfo();

using namespace amrex;
using namespace std::chrono;

void add_par () {
   ParmParse pp("eb2");
   if(not pp.contains("extend_domain_face")) {
      pp.add("extend_domain_face",true);
   }
}

int main(int argc, char* argv[])
{

    // check to see if it contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                writeBuildInfo();
                return 0;
            }
        }
    }

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, add_par);
    { /* These braces are necessary to ensure amrex::Finalize() can be called without explicitly
        deleting all the incflo member MultiFabs */

        BL_PROFILE("main()");

        // Issue an error if input file is not given
        if(argc < 2) amrex::Abort("Input file must be given as command-line argument.");

        // Write out the incflo git hash (the AMReX git hash is already written)
        const char* githash_incflo = buildInfoGetGitHash(1);
        amrex::Print() << "incflo git hash: " << githash_incflo << "\n";

        // Start timing the program
        Real start_time = Real(ParallelDescriptor::second());

        // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh.
        incflo my_incflo;

        // Initialize data, parameters, arrays and derived internals
        my_incflo.InitData();
        
        if (my_incflo.is_restart_run()) {

            if (my_incflo.get_seed() > 0) {
                // initializes the seed for C++ random number calls
                InitRandom(my_incflo.get_seed()+ParallelDescriptor::MyProc(),
                           ParallelDescriptor::NProcs(),
                           my_incflo.get_seed()+ParallelDescriptor::MyProc());
            } else if (my_incflo.get_seed() == 0) {
                // initializes the seed for C++ random number calls based on the clock
                auto now = time_point_cast<nanoseconds>(system_clock::now());
                int randSeed = now.time_since_epoch().count();
                // broadcast the same root seed to all processors
                ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
                InitRandom(randSeed+ParallelDescriptor::MyProc(),
                           ParallelDescriptor::NProcs(),
                           randSeed+ParallelDescriptor::MyProc());
            } else {
                Abort("Must supply non-negative seed");
            }
        }

        // Time spent on initialization
        Real init_time = Real(ParallelDescriptor::second()) - start_time;

        // Evolve system to final time
        my_incflo.Evolve();

        // Time spent in total
        Real end_time = Real(ParallelDescriptor::second()) - start_time;

        ParallelDescriptor::ReduceRealMax(init_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

        // Print timing results
        amrex::Print() << "Time spent in InitData():    " << init_time << "\n";
        amrex::Print() << "Time spent in Evolve():      " << end_time - init_time << "\n";
    }
    amrex::Finalize();
}
