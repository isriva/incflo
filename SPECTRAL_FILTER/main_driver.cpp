#include "spectral_functions.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_String.H>

using namespace amrex;

namespace {

std::string PreviousCheckpointName(const std::string& restart_file, int step)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        step > 0,
        "restart_file must identify a checkpoint after step zero");

    std::string const step_suffix = amrex::Concatenate("", step);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        restart_file.size() >= step_suffix.size() &&
            restart_file.compare(restart_file.size() - step_suffix.size(),
                                 step_suffix.size(), step_suffix) == 0,
        "restart_file name does not match its checkpoint Header step");

    std::string const checkpoint_root =
        restart_file.substr(0, restart_file.size() - step_suffix.size());
    return amrex::Concatenate(checkpoint_root, step - 1);
}

void AssertMatchingCheckpointGeometry(const Geometry& geom,
                                      const BoxArray& ba,
                                      const Geometry& previous_geom,
                                      const BoxArray& previous_ba)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        previous_ba == ba,
        "Previous checkpoint grid does not match the current checkpoint grid");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        previous_geom.Domain() == geom.Domain(),
        "Previous checkpoint domain does not match the current checkpoint domain");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        previous_geom.Coord() == geom.Coord(),
        "Previous checkpoint coordinate system does not match the current checkpoint");
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            previous_geom.ProbLo(idim) == geom.ProbLo(idim) &&
                previous_geom.ProbHi(idim) == geom.ProbHi(idim) &&
                previous_geom.isPeriodic(idim) == geom.isPeriodic(idim),
            "Previous checkpoint physical geometry does not match the current checkpoint");
    }
}

} // namespace

void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()", main_driver);
    amrex::ignore_unused(argv);

    Real ts1 = ParallelDescriptor::second();

    ParmParse pp;

    std::string restart_file;
    if (!pp.query("restart_file", restart_file)) {
        Abort("SPECTRAL_FILTER requires restart_file=<incflo checkpoint directory>");
    }

    Real kmin = 0.0;
    // Real kmax = 0.0;
    std::vector<Real> kmax_list;
    if (!pp.query("kmin", kmin)) {
        Abort("SPECTRAL_FILTER requires kmin=<minimum integer wavenumber>");
    }
    // if (!pp.query("kmax", kmax)) {
    //     Abort("SPECTRAL_FILTER requires kmax=<maximum integer wavenumber>");
    // }
    if (!pp.queryarr("kmax_list", kmax_list)) {
        Abort("SPECTRAL_FILTER requires kmax_list=<list of maximum integer wavenumbers>");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmin >= 0.0, "kmin must be non-negative");
    // AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmax >= kmin, "kmax must be greater than or equal to kmin");
    // Loop through the vector to validate each kmax value
    for (Real kmax : kmax_list) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmax >= kmin, "Every kmax must be greater than or equal to kmin");
    }

    std::string filter_type_name = "sharp";
    pp.query("filter_type", filter_type_name);
    SpectralFilterOptions filter_options;
    if (filter_type_name == "sharp") {
        filter_options.filter_type = SpectralFilterType::Sharp;
    } else if (filter_type_name == "sinc_sq") {
        filter_options.filter_type = SpectralFilterType::SincSq;
    } else {
        Abort("filter_type must be either sharp or sinc_sq");
    }

    int zero_outside_range = 1;
    pp.query("zero_outside_range", zero_outside_range);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        zero_outside_range == 0 || zero_outside_range == 1,
        "zero_outside_range must be either 0 or 1");
    filter_options.zero_outside_range = (zero_outside_range == 1);

    if (filter_options.filter_type == SpectralFilterType::SincSq) {
        for (Real kmax : kmax_list) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmax > 0.0, "sinc_sq requires every kmax to be positive");
        }
    }
    int plot_filter = 0;
    pp.query("plot_filter", plot_filter);

    int plot_fourier = 0;
    pp.query("plot_fourier", plot_fourier);

    int use_prime_tau_int = 0;
    pp.query("use_prime_tau", use_prime_tau_int);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        use_prime_tau_int == 0 || use_prime_tau_int == 1,
        "use_prime_tau must be either 0 or 1");
    bool const use_prime_tau = (use_prime_tau_int == 1);

    Vector<int> is_periodic(AMREX_SPACEDIM, 1);
    int const nper = pp.countval("is_periodic");
    if (nper > 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nper == AMREX_SPACEDIM,
            "is_periodic must contain AMREX_SPACEDIM entries");
        pp.getarr("is_periodic", is_periodic, 0, AMREX_SPACEDIM);
    }

    Geometry geom;
    BoxArray ba;
    DistributionMapping dmap;
    MultiFab velocity;
    int step = 0;
    Real time = 0.0;

    SpectralReadCheckPoint(restart_file, is_periodic, geom, ba, dmap, velocity, step, time);

    std::string const previous_restart_file = PreviousCheckpointName(restart_file, step);
    Geometry previous_geom;
    BoxArray previous_ba;
    DistributionMapping previous_dmap;
    MultiFab previous_velocity;
    int previous_step = 0;
    Real previous_time = 0.0;
    SpectralReadCheckPoint(previous_restart_file, is_periodic, previous_geom, previous_ba,
                           previous_dmap, previous_velocity, previous_step, previous_time);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        previous_step == step - 1,
        "Previous checkpoint Header step must equal the current Header step minus one");
    AssertMatchingCheckpointGeometry(geom, ba, previous_geom, previous_ba);
    amrex::ignore_unused(previous_time);

    MultiFab velocity_filter(ba, dmap, 3, 0);
    MultiFab previous_velocity_filter(previous_ba, previous_dmap, 3, 0);

#if (AMREX_SPACEDIM == 3)
    int constexpr num_vv_comps = 6;
#else
    int constexpr num_vv_comps = 3;
#endif

    MultiFab vv_filter(ba, dmap, num_vv_comps, 0);
    

    velocity.FillBoundary(geom.periodicity());
    previous_velocity.FillBoundary(previous_geom.periodicity());

    // Loop over all of the kmax values
    for (Real kmax : kmax_list) {
        velocity_filter.setVal(0.0);
        previous_velocity_filter.setVal(0.0);
        vv_filter.setVal(0.0);
        
        // Filter both checkpoints with identical spectral bounds.
        SpectralVelDecomp(velocity, velocity_filter, kmin, kmax, filter_options, geom);
        SpectralVelDecomp(previous_velocity, previous_velocity_filter, kmin, kmax,
                          filter_options, previous_geom);
        velocity_filter.FillBoundary(geom.periodicity());
        previous_velocity_filter.FillBoundary(previous_geom.periodicity());

        // Filter the outer product of the velocity
        if (use_prime_tau) {
            MultiFab velocity_prime(ba, dmap, AMREX_SPACEDIM, 0);
            MultiFab::Copy(velocity_prime, velocity, 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Subtract(velocity_prime, velocity_filter, 0, 0, AMREX_SPACEDIM, 0);
            SpectralVelProductDecomp(velocity_prime, vv_filter, kmin, kmax,
                                     filter_options, geom);
        } else {
            SpectralVelProductDecomp(velocity, vv_filter, kmin, kmax,
                                     filter_options, geom);
        }
        vv_filter.FillBoundary(geom.periodicity());
    
        // velocity.FillBoundary(geom.periodicity());
        // SpectralVelDecomp(velocity, velocity_filter, kmin, kmax, geom);
        // velocity_filter.FillBoundary(geom.periodicity());
    
        if (plot_filter != 0) {
            SpectralWritePlotFile(
                step, kmin, kmax, filter_options, geom, velocity, velocity_filter,
                previous_velocity_filter, vv_filter, use_prime_tau);
        }
        if (plot_fourier != 0) {
            SpectralWriteFourierPlotFile(step, kmin, kmax, filter_options, geom,
                                         velocity_filter, vv_filter);
        }
    }

    Real ts2 = ParallelDescriptor::second() - ts1;
    ParallelDescriptor::ReduceRealMax(ts2, ParallelDescriptor::IOProcessorNumber());
    Print() << "Time (spectral filtering) " << ts2 << " seconds\n";

    Long min_fab_megabytes = TotalBytesAllocatedInFabsHWM() / 1048576;
    Long max_fab_megabytes = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, ParallelDescriptor::IOProcessorNumber());

    Print() << "High-water FAB megabyte spread across MPI nodes: ["
            << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
}
