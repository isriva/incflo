#include "spectral_functions.H"

#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <iomanip>
#include <sstream>

namespace {

void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void ReadProbLine(std::istream& is, amrex::Real* prob, const char* name)
{
    std::string line;
    std::getline(is, line);
    std::istringstream lis(line);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (!(lis >> prob[idim])) {
            amrex::Abort(std::string("Checkpoint Header is missing ") + name + " entries");
        }
    }
}

std::string FilterPlotFileName(int step, amrex::Real kmin, amrex::Real kmax)
{
    std::ostringstream os;
    os << "filtered_" << step << "_" << std::setprecision(12) << kmin << "_" << kmax;
    return os.str();
}

std::string FourierPlotFileName(int step, amrex::Real kmin, amrex::Real kmax)
{
    std::ostringstream os;
    os << "filtered_fourier_" << step << "_" << std::setprecision(12) << kmin << "_" << kmax;
    return os.str();
}

void FillVelocityWithGhosts(const amrex::MultiFab& velocity,
                            amrex::MultiFab& velocity_g,
                            const amrex::Geometry& geom)
{
    velocity_g.setVal(0.0);
    amrex::MultiFab::Copy(velocity_g, velocity, 0, 0, velocity.nComp(), 0);
    velocity_g.FillBoundary(geom.periodicity());
}

void ComputeVorticityField(const amrex::MultiFab& velocity_g,
                           amrex::MultiFab& vort,
                           const amrex::Geometry& geom)
{
    amrex::Real const idx = geom.InvCellSize(0);
    amrex::Real const idy = geom.InvCellSize(1);
#if (AMREX_SPACEDIM == 3)
    amrex::Real const idz = geom.InvCellSize(2);
#endif

    for (amrex::MFIter mfi(vort, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& vel = velocity_g.const_array(mfi);
        amrex::Array4<amrex::Real> const& vo = vort.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real const vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
            amrex::Real const uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
#if (AMREX_SPACEDIM == 2)
            vo(i,j,k) = vx - uy;
#else
            amrex::Real const wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;
            amrex::Real const wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;
            amrex::Real const uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
            amrex::Real const vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
            vo(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }
}

void ComputeFullDomainFFT(const amrex::MultiFab& variables,
                          amrex::MultiFab& variables_dft_real,
                          amrex::MultiFab& variables_dft_imag)
{
    amrex::Box const domain = variables.boxArray().minimalBox();

    bool const chopped_in_x = domain.length(0) > 1;
    bool const chopped_in_y = !chopped_in_x && domain.length(1) > 1;
#if (AMREX_SPACEDIM == 3)
    bool const chopped_in_z = !chopped_in_x && !chopped_in_y && domain.length(2) > 1;
#endif
#if (AMREX_SPACEDIM == 3)
    bool const can_transform = chopped_in_x || chopped_in_y || chopped_in_z;
#else
    bool const can_transform = chopped_in_x || chopped_in_y;
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        can_transform,
        "ComputeFullDomainFFT: cannot transform a one-cell domain");

    amrex::Real const sqrtnpts = std::sqrt(amrex::Real(domain.numPts()));

    amrex::BoxArray const ba = variables.boxArray();
    amrex::DistributionMapping const dm = variables.DistributionMap();
    amrex::MultiFab phi(ba, dm, 1, 0);

    amrex::BoxArray ba_onegrid(domain);
    amrex::DistributionMapping dm_onegrid(ba_onegrid);

    amrex::FFT::R2C<amrex::Real, amrex::FFT::Direction::forward> fft(domain);
    auto const& [ba_fft, dm_fft] = fft.getSpectralDataLayout();
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<amrex::Real>>> phi_fft(ba_fft, dm_fft, 1, 0);

    amrex::Box const domain_fft = ba_fft.minimalBox();
    amrex::BoxArray ba_fft_onegrid(domain_fft);
    amrex::FabArray<amrex::BaseFab<amrex::GpuComplex<amrex::Real>>> phi_fft_onegrid(ba_fft_onegrid, dm_onegrid, 1, 0);

    amrex::MultiFab real_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    amrex::MultiFab imag_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    for (int comp = 0; comp < variables.nComp(); ++comp) {
        amrex::MultiFab::Copy(phi, variables, comp, 0, 1, 0);
        fft.forward(phi, phi_fft);
        phi_fft_onegrid.ParallelCopy(phi_fft, 0, 0, 1);

        for (amrex::MFIter mfi(real_onegrid); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.fabbox();
            amrex::Array4<const amrex::GpuComplex<amrex::Real>> const& spectral = phi_fft_onegrid.const_array(mfi);
            amrex::Array4<amrex::Real> const& realpart = real_onegrid.array(mfi);
            amrex::Array4<amrex::Real> const& imagpart = imag_onegrid.array(mfi);

            if (chopped_in_x) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i <= bx.length(0) / 2) {
                        realpart(i,j,k) = spectral(i,j,k).real();
                        imagpart(i,j,k) = spectral(i,j,k).imag();
                    } else {
                        int const iloc = bx.length(0) - i;
                        int const jloc = (j == 0) ? 0 : bx.length(1) - j;
#if (AMREX_SPACEDIM == 2)
                        int const kloc = 0;
#else
                        int const kloc = (k == 0) ? 0 : bx.length(2) - k;
#endif
                        realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
                        imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
            } else if (chopped_in_y) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j <= bx.length(1) / 2) {
                        realpart(i,j,k) = spectral(i,j,k).real();
                        imagpart(i,j,k) = spectral(i,j,k).imag();
                    } else {
                        int const iloc = (i == 0) ? 0 : bx.length(0) - i;
                        int const jloc = bx.length(1) - j;
#if (AMREX_SPACEDIM == 2)
                        int const kloc = 0;
#else
                        int const kloc = (k == 0) ? 0 : bx.length(2) - k;
#endif
                        realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
                        imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
#if (AMREX_SPACEDIM == 3)
            } else if (chopped_in_z) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k <= bx.length(2) / 2) {
                        realpart(i,j,k) = spectral(i,j,k).real();
                        imagpart(i,j,k) = spectral(i,j,k).imag();
                    } else {
                        int const iloc = (i == 0) ? 0 : bx.length(0) - i;
                        int const jloc = (j == 0) ? 0 : bx.length(1) - j;
                        int const kloc = bx.length(2) - k;
                        realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
                        imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
#endif
            }
        }

        variables_dft_real.ParallelCopy(real_onegrid, 0, comp, 1);
        variables_dft_imag.ParallelCopy(imag_onegrid, 0, comp, 1);
    }
}

} // namespace

void SpectralReadCheckPoint(const std::string& restart_file,
                            const amrex::Vector<int>& is_periodic,
                            amrex::Geometry& geom,
                            amrex::BoxArray& ba,
                            amrex::DistributionMapping& dmap,
                            amrex::MultiFab& velocity,
                            int& step,
                            amrex::Real& time)
{
    BL_PROFILE_VAR("SpectralReadCheckPoint()", SpectralReadCheckPoint);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        is_periodic.size() == AMREX_SPACEDIM,
        "is_periodic must contain AMREX_SPACEDIM entries");

    amrex::Print() << "Restarting from checkpoint " << restart_file << "\n";

    std::string const header_file = restart_file + "/Header";
    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(header_file, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line;

    std::getline(is, line); // title

    int finest_level = -1;
    is >> finest_level;
    GotoNextLine(is);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        finest_level == 0,
        "SPECTRAL_FILTER supports only single-level incflo checkpoints (finest_level must be 0)");

    is >> step;
    GotoNextLine(is);

    is >> time;
    GotoNextLine(is);

    amrex::Real dt = 0.0;
    amrex::Real prev_dt = 0.0;
    amrex::Real prev_prev_dt = 0.0;
    is >> dt;
    GotoNextLine(is);
    is >> prev_dt;
    GotoNextLine(is);
    is >> prev_prev_dt;
    GotoNextLine(is);
    amrex::ignore_unused(dt, prev_dt, prev_prev_dt);

    amrex::Real prob_lo[AMREX_SPACEDIM];
    amrex::Real prob_hi[AMREX_SPACEDIM];
    ReadProbLine(is, prob_lo, "prob_lo");
    ReadProbLine(is, prob_hi, "prob_hi");

    ba.readFrom(is);
    GotoNextLine(is);

    amrex::Box const domain = ba.minimalBox();
    amrex::RealBox real_box(prob_lo, prob_hi);
    geom.define(domain, &real_box, amrex::CoordSys::cartesian, is_periodic.data());

    dmap.define(ba, amrex::ParallelDescriptor::NProcs());

    amrex::VisMF::Read(
        velocity,
        amrex::MultiFabFileFullPrefix(0, restart_file, "Level_", "velocity"));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity.nComp() >= AMREX_SPACEDIM,
        "Checkpoint velocity has fewer than AMREX_SPACEDIM components");

    ba = velocity.boxArray();
    dmap = velocity.DistributionMap();

    amrex::Print() << "Read level 0 velocity with " << velocity.nComp()
                   << " components on domain " << geom.Domain() << "\n";
}

void SpectralVelDecomp(const amrex::MultiFab& velocity,
                       amrex::MultiFab& velocity_filter,
                       amrex::Real kmin,
                       amrex::Real kmax,
                       const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("SpectralVelDecomp()", SpectralVelDecomp);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity.nComp() >= AMREX_SPACEDIM,
        "SpectralVelDecomp: input velocity must have at least AMREX_SPACEDIM components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity_filter.nComp() == 3,
        "SpectralVelDecomp: velocity_filter must have exactly 3 components");

    amrex::Box const domain = geom.Domain();
    int const nx = domain.length(0);
#if (AMREX_SPACEDIM >= 2)
    int const ny = domain.length(1);
#else
    int const ny = 1;
#endif
#if (AMREX_SPACEDIM == 3)
    int const nz = domain.length(2);
#else
    int const nz = 1;
#endif

    amrex::Real const kmin2 = kmin * kmin;
    amrex::Real const kmax2 = kmax * kmax;

    amrex::FFT::R2C<> fft(domain);
    amrex::Real const scale = fft.scalingFactor();

    amrex::MultiFab velocity_single(velocity.boxArray(), velocity.DistributionMap(), 1, 0);
    amrex::MultiFab filtered_single(velocity.boxArray(), velocity.DistributionMap(), 1, 0);

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        amrex::MultiFab::Copy(velocity_single, velocity, comp, 0, 1, 0);
        filtered_single.setVal(0.0);

        fft.forwardThenBackward(
            velocity_single,
            filtered_single,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::GpuComplex<amrex::Real>& sp) noexcept
            {
                int const ik = (i <= nx / 2) ? i : i - nx;
                int const jk = (j <= ny / 2) ? j : j - ny;
#if (AMREX_SPACEDIM == 3)
                int const kk = (k <= nz / 2) ? k : k - nz;
#else
                amrex::ignore_unused(k, nz);
                int const kk = 0;
#endif
                amrex::Real const ksq = amrex::Real(ik*ik + jk*jk + kk*kk);
                if (ksq < kmin2 || ksq > kmax2) {
                    sp = amrex::GpuComplex<amrex::Real>(0.0, 0.0);
                } else {
                    sp *= scale;
                }
            });

        amrex::MultiFab::Copy(velocity_filter, filtered_single, 0, comp, 1, 0);
    }

#if (AMREX_SPACEDIM < 3)
    velocity_filter.setVal(0.0, 2, 1, 0);
#endif
}

void SpectralWritePlotFile(int step,
                           amrex::Real kmin,
                           amrex::Real kmax,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& velocity,
                           const amrex::MultiFab& velocity_filter)
{
    BL_PROFILE_VAR("SpectralWritePlotFile()", SpectralWritePlotFile);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity.nComp() >= AMREX_SPACEDIM,
        "SpectralWritePlotFile: input velocity must have at least AMREX_SPACEDIM components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity_filter.nComp() == 3,
        "SpectralWritePlotFile: velocity_filter must have exactly 3 components");

#if (AMREX_SPACEDIM == 2)
    int constexpr nplot = 6;
#else
    int constexpr nplot = 8;
#endif

    amrex::MultiFab velocity_g(velocity.boxArray(), velocity.DistributionMap(), 3, 1);
    amrex::MultiFab velocity_filter_g(velocity_filter.boxArray(), velocity_filter.DistributionMap(), 3, 1);
    velocity_g.setVal(0.0);
    velocity_filter_g.setVal(0.0);
    amrex::MultiFab::Copy(velocity_g, velocity, 0, 0, AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(velocity_filter_g, velocity_filter, 0, 0, 3, 0);
    velocity_g.FillBoundary(geom.periodicity());
    velocity_filter_g.FillBoundary(geom.periodicity());

    amrex::MultiFab output(velocity.boxArray(), velocity.DistributionMap(), nplot, 0);
    output.setVal(0.0);

    amrex::Real const idx = geom.InvCellSize(0);
    amrex::Real const idy = geom.InvCellSize(1);
#if (AMREX_SPACEDIM == 3)
    amrex::Real const idz = geom.InvCellSize(2);
#endif

    for (amrex::MFIter mfi(output, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& out = output.array(mfi);
        amrex::Array4<const amrex::Real> const& vel = velocity_g.const_array(mfi);
        amrex::Array4<const amrex::Real> const& filt = velocity_filter_g.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            out(i,j,k,0) = vel(i,j,k,0);
            out(i,j,k,1) = vel(i,j,k,1);
#if (AMREX_SPACEDIM == 2)
            out(i,j,k,2) = filt(i,j,k,0);
            out(i,j,k,3) = filt(i,j,k,1);

            amrex::Real const vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
            amrex::Real const uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            amrex::Real const vx_f = amrex::Real(0.5) * (filt(i+1,j,k,1) - filt(i-1,j,k,1)) * idx;
            amrex::Real const uy_f = amrex::Real(0.5) * (filt(i,j+1,k,0) - filt(i,j-1,k,0)) * idy;
            out(i,j,k,4) = vx - uy;
            out(i,j,k,5) = vx_f - uy_f;
#else
            out(i,j,k,2) = vel(i,j,k,2);
            out(i,j,k,3) = filt(i,j,k,0);
            out(i,j,k,4) = filt(i,j,k,1);
            out(i,j,k,5) = filt(i,j,k,2);

            amrex::Real const vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
            amrex::Real const wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;
            amrex::Real const uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            amrex::Real const wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;
            amrex::Real const uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
            amrex::Real const vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;

            amrex::Real const vx_f = amrex::Real(0.5) * (filt(i+1,j,k,1) - filt(i-1,j,k,1)) * idx;
            amrex::Real const wx_f = amrex::Real(0.5) * (filt(i+1,j,k,2) - filt(i-1,j,k,2)) * idx;
            amrex::Real const uy_f = amrex::Real(0.5) * (filt(i,j+1,k,0) - filt(i,j-1,k,0)) * idy;
            amrex::Real const wy_f = amrex::Real(0.5) * (filt(i,j+1,k,2) - filt(i,j-1,k,2)) * idy;
            amrex::Real const uz_f = amrex::Real(0.5) * (filt(i,j,k+1,0) - filt(i,j,k-1,0)) * idz;
            amrex::Real const vz_f = amrex::Real(0.5) * (filt(i,j,k+1,1) - filt(i,j,k-1,1)) * idz;

            out(i,j,k,6) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
            out(i,j,k,7) = std::sqrt((wy_f-vz_f)*(wy_f-vz_f) +
                                      (uz_f-wx_f)*(uz_f-wx_f) +
                                      (vx_f-uy_f)*(vx_f-uy_f));
#endif
        });
    }

#if (AMREX_SPACEDIM == 2)
    amrex::Vector<std::string> varNames{
        "velx", "vely",
        "velx_filter", "vely_filter",
        "vort", "vort_filter"};
#else
    amrex::Vector<std::string> varNames{
        "velx", "vely", "velz",
        "velx_filter", "vely_filter", "velz_filter",
        "vort", "vort_filter"};
#endif

    std::string const plotfilename = FilterPlotFileName(step, kmin, kmax);
    amrex::Print() << "Writing filtered velocity plotfile " << plotfilename << "\n";
    amrex::WriteSingleLevelPlotfile(plotfilename, output, varNames, geom, 0.0, step);
}


void SpectralWriteFourierPlotFile(int step,
                                  amrex::Real kmin,
                                  amrex::Real kmax,
                                  const amrex::Geometry& geom,
                                  const amrex::MultiFab& velocity_filter)
{
    BL_PROFILE_VAR("SpectralWriteFourierPlotFile()", SpectralWriteFourierPlotFile);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity_filter.nComp() == 3,
        "SpectralWriteFourierPlotFile: velocity_filter must have exactly 3 components");

    int constexpr nfields = AMREX_SPACEDIM + 1;

    amrex::MultiFab velocity_filter_g(velocity_filter.boxArray(), velocity_filter.DistributionMap(), 3, 1);
    FillVelocityWithGhosts(velocity_filter, velocity_filter_g, geom);

    amrex::MultiFab fields(velocity_filter.boxArray(), velocity_filter.DistributionMap(), nfields, 0);
    fields.setVal(0.0);
    amrex::MultiFab::Copy(fields, velocity_filter, 0, 0, AMREX_SPACEDIM, 0);

    amrex::MultiFab vort_filter(fields, amrex::make_alias, AMREX_SPACEDIM, 1);
    ComputeVorticityField(velocity_filter_g, vort_filter, geom);

    amrex::MultiFab fourier_real(fields.boxArray(), fields.DistributionMap(), nfields, 0);
    amrex::MultiFab fourier_imag(fields.boxArray(), fields.DistributionMap(), nfields, 0);
    fourier_real.setVal(0.0);
    fourier_imag.setVal(0.0);
    ComputeFullDomainFFT(fields, fourier_real, fourier_imag);

    amrex::MultiFab output(fields.boxArray(), fields.DistributionMap(), 2*nfields, 0);
    output.setVal(0.0);
    amrex::MultiFab::Copy(output, fourier_real, 0, 0, nfields, 0);
    amrex::MultiFab::Copy(output, fourier_imag, 0, nfields, nfields, 0);

#if (AMREX_SPACEDIM == 2)
    amrex::Vector<std::string> baseNames{
        "velx_filter", "vely_filter", "vort_filter"};
#else
    amrex::Vector<std::string> baseNames{
        "velx_filter", "vely_filter", "velz_filter", "vort_filter"};
#endif

    amrex::Vector<std::string> varNames(2*nfields);
    for (int n = 0; n < nfields; ++n) {
        varNames[n] = baseNames[n] + "_real";
        varNames[n+nfields] = baseNames[n] + "_imag";
    }

    std::string const plotfilename = FourierPlotFileName(step, kmin, kmax);
    amrex::Print() << "Writing filtered Fourier plotfile " << plotfilename << "\n";
    amrex::WriteSingleLevelPlotfile(plotfilename, output, varNames, geom, 0.0, step);
}
