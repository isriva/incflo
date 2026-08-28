#include "spectral_functions.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <cmath>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr int N = 32;
        constexpr int mode = 10;
        constexpr amrex::Real two_pi = 2.0 * amrex::Math::pi<amrex::Real>();

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(N - 1));
        amrex::BoxArray ba(domain);
        amrex::DistributionMapping dmap(ba);
        amrex::RealBox real_box({0.0, 0.0}, {two_pi, two_pi});
        amrex::Array<int, AMREX_SPACEDIM> periodic{{1, 1}};
        amrex::Geometry geom(domain, real_box, 0, periodic);

        amrex::MultiFab velocity(ba, dmap, 2, 0);
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.tilebox();
            auto const& vel = velocity.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::Real const u = std::sin(two_pi * mode * i / N);
                vel(i,j,k,0) = u;
                vel(i,j,k,1) = u;
            });
        }

        amrex::MultiFab vv_filter(ba, dmap, 3, 0);
        // 2*mode=20 is above the retained cutoff and aliases to 12 without
        // dealiasing. The correct filtered result is sin(mode*x)^2 -> 1/2.
        SpectralVelProductDecomp(velocity, vv_filter, 0.0, 15.0, geom);

        amrex::Real max_error = amrex::ReduceMax(
            vv_filter, 0,
            [] AMREX_GPU_HOST_DEVICE (amrex::Box const& bx,
                                       amrex::Array4<const amrex::Real> const& vv) noexcept
            {
                amrex::Real error = 0.0;
                amrex::Loop(bx, [&] (int i, int j, int k) noexcept
                {
                    error = amrex::max(error, amrex::Math::abs(vv(i,j,k,0) - 0.5));
                });
                return error;
            });

        amrex::ParallelDescriptor::ReduceRealMax(
            max_error, amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "max error = " << max_error << "\n";
            if (max_error > 1.0e-10) {
                amrex::Abort("SpectralVelProductDecomp dealiasing test failed");
            }
            amrex::Print() << "PASS\n";
        }
    }
    amrex::Finalize();
    return 0;
}
