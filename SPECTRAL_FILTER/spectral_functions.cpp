#include "spectral_functions.H"

#include "IncfloStructFact.H"

#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <set>
#include <string>

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

std::string Delta_nu_FileName(int step)
{
    return amrex::Concatenate("Delta_nu_", step, 7) + ".txt";
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

using ComplexFab = amrex::BaseFab<amrex::GpuComplex<amrex::Real>>;
using ComplexFabArray = amrex::FabArray<ComplexFab>;

// Copy Fourier modes between real-to-complex spectra with different real-space
// sizes. The first spectral direction is the nonnegative half-spectrum; the
// remaining directions contain both positive and negative modes.
void RemapFourierModes(const ComplexFabArray& source,
                       ComplexFabArray& destination,
                       const amrex::IntVect& source_size,
                       const amrex::IntVect& destination_size)
{
    destination.setVal(amrex::GpuComplex<amrex::Real>(0.0, 0.0));

    amrex::MFIter dmfi(destination);
    amrex::MFIter smfi(source);
    for (; dmfi.isValid(); ++dmfi, ++smfi) {
        amrex::Box const& dbx = dmfi.fabbox();
        amrex::Box const& sbx = smfi.fabbox();
        auto const& src = source.const_array(smfi);
        auto const& dst = destination.array(dmfi);

        amrex::ParallelFor(dbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int const dcoord[3] = {i, j, k};
            int scoord[3] = {0, 0, 0};
            bool valid = true;

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                int const dmode = dcoord[idim] - dbx.smallEnd(idim);
                if (idim == 0) {
                    valid = dmode <= source_size[idim] / 2;
                    scoord[idim] = sbx.smallEnd(idim) + dmode;
                } else {
                    int mode = dmode;
                    if (mode > destination_size[idim] / 2) {
                        mode -= destination_size[idim];
                    }
                    int const abs_mode = (mode < 0) ? -mode : mode;
                    valid = valid && (abs_mode <= source_size[idim] / 2);
                    int const source_mode = (mode < 0) ? mode + source_size[idim] : mode;
                    scoord[idim] = sbx.smallEnd(idim) + source_mode;
                }
            }

            if (valid) {
                dst(i,j,k) = src(scoord[0], scoord[1], scoord[2]);
            }
        });
    }
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



void SpectralVelProductDecomp(const amrex::MultiFab& velocity,
                              amrex::MultiFab& vv_filter,
                              amrex::Real kmin,
                              amrex::Real kmax,
                              const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("SpectralVelProductDecomp()", SpectralVelProductDecomp);

#if (AMREX_SPACEDIM == 3)
    int constexpr num_vv_comps = 6; // 0:uu, 1:vv, 2:ww, 3:uv, 4:uw, 5:vw
#else
    int constexpr num_vv_comps = 3; // 0:uu, 1:vv, 2:uv
#endif

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity.nComp() >= AMREX_SPACEDIM,
        "SpectralVelProductDecomp: input velocity must have at least AMREX_SPACEDIM components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        vv_filter.nComp() >= num_vv_comps,
        "SpectralVelProductDecomp: vv_filter must have enough components to store the symmetric tensor");

    constexpr amrex::Real padding_factor = amrex::Real(1.5);
    amrex::Box const domain = geom.Domain();
    amrex::IntVect original_size = amrex::IntVect::TheZeroVector();
    amrex::IntVect padded_size = amrex::IntVect::TheZeroVector();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const n = domain.length(idim);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            n % 2 == 0,
            "SpectralVelProductDecomp: every domain dimension must be even for 3/2 padding");
        original_size[idim] = n;
        padded_size[idim] = static_cast<int>(padding_factor * n);
    }

    amrex::Box padded_domain = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        padded_domain.setBig(idim, padded_domain.smallEnd(idim) + padded_size[idim] - 1);
    }

    amrex::BoxArray const ba_onegrid(domain);
    amrex::DistributionMapping const dm_onegrid(ba_onegrid);
    amrex::BoxArray const ba_padded_onegrid(padded_domain);
    amrex::DistributionMapping const dm_padded_onegrid(ba_padded_onegrid);

    amrex::MultiFab velocity_onegrid(ba_onegrid, dm_onegrid, AMREX_SPACEDIM, 0);
    velocity_onegrid.ParallelCopy(velocity, 0, 0, AMREX_SPACEDIM, 0, 0);

    amrex::FFT::R2C<> original_fft(domain);
    amrex::FFT::R2C<> padded_fft(padded_domain);
    auto const& [ba_original_fft, dm_original_fft] = original_fft.getSpectralDataLayout();
    auto const& [ba_padded_fft, dm_padded_fft] = padded_fft.getSpectralDataLayout();

    amrex::FabArray<ComplexFab> original_fft_dist(ba_original_fft, dm_original_fft, 1, 0);
    amrex::FabArray<ComplexFab> padded_fft_dist(ba_padded_fft, dm_padded_fft, 1, 0);
    amrex::FabArray<ComplexFab> original_fft_onegrid(
        amrex::BoxArray(ba_original_fft.minimalBox()), dm_onegrid, 1, 0);
    amrex::FabArray<ComplexFab> padded_fft_onegrid(
        amrex::BoxArray(ba_padded_fft.minimalBox()), dm_padded_onegrid, 1, 0);

    amrex::MultiFab padded_product(ba_padded_onegrid, dm_padded_onegrid, num_vv_comps, 0);
    amrex::MultiFab filtered_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    amrex::Real const padded_scale = padded_fft.scalingFactor();
    amrex::Real const original_scale = original_fft.scalingFactor();
    amrex::Real const kmin2 = kmin * kmin;
    amrex::Real const kmax2 = kmax * kmax;

    // Dealias the velocity by padding its spectrum, then transform once to the
    // padded real-space grid.
    amrex::MultiFab padded_velocity_all(ba_padded_onegrid, dm_padded_onegrid,
                                         AMREX_SPACEDIM, 0);
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        amrex::MultiFab velocity_component(ba_onegrid, dm_onegrid, 1, 0);
        amrex::MultiFab::Copy(velocity_component, velocity_onegrid, comp, 0, 1, 0);
        original_fft.forward(velocity_component, original_fft_dist);
        original_fft_onegrid.ParallelCopy(original_fft_dist, 0, 0, 1);
        RemapFourierModes(original_fft_onegrid, padded_fft_onegrid,
                          original_size, padded_size);
        padded_fft_dist.ParallelCopy(padded_fft_onegrid, 0, 0, 1);
        for (amrex::MFIter mfi(padded_fft_dist); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.fabbox();
            auto const& spectrum = padded_fft_dist.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                spectrum(i,j,k) *= padded_scale;
            });
        }
        amrex::MultiFab component_padded(ba_padded_onegrid, dm_padded_onegrid, 1, 0);
        padded_fft.backward(padded_fft_dist, component_padded);
        amrex::MultiFab::Copy(padded_velocity_all, component_padded, 0, comp, 1, 0);
    }

    for (amrex::MFIter mfi(padded_product, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& bx = mfi.tilebox();
        auto const& vel = padded_velocity_all.const_array(mfi);
        auto const& product = padded_product.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            product(i,j,k,0) = vel(i,j,k,0) * vel(i,j,k,0);
            product(i,j,k,1) = vel(i,j,k,1) * vel(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
            product(i,j,k,2) = vel(i,j,k,2) * vel(i,j,k,2);
            product(i,j,k,3) = vel(i,j,k,0) * vel(i,j,k,1);
            product(i,j,k,4) = vel(i,j,k,0) * vel(i,j,k,2);
            product(i,j,k,5) = vel(i,j,k,1) * vel(i,j,k,2);
#else
            product(i,j,k,2) = vel(i,j,k,0) * vel(i,j,k,1);
#endif
        });
    }

    for (int comp = 0; comp < num_vv_comps; ++comp) {
        amrex::MultiFab product_component(ba_padded_onegrid, dm_padded_onegrid, 1, 0);
        amrex::MultiFab::Copy(product_component, padded_product, comp, 0, 1, 0);
        padded_fft.forward(product_component, padded_fft_dist);
        padded_fft_onegrid.ParallelCopy(padded_fft_dist, 0, 0, 1);

        for (amrex::MFIter mfi(padded_fft_onegrid); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.fabbox();
            auto const& spectrum = padded_fft_onegrid.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int const ik = (i <= padded_size[0] / 2) ? i : i - padded_size[0];
                int const jk = (j <= padded_size[1] / 2) ? j : j - padded_size[1];
#if (AMREX_SPACEDIM == 3)
                int const kk = (k <= padded_size[2] / 2) ? k : k - padded_size[2];
#else
                amrex::ignore_unused(k);
                int const kk = 0;
#endif
                amrex::Real const ksq = amrex::Real(ik*ik + jk*jk + kk*kk);
                if (ksq < kmin2 || ksq > kmax2) {
                    spectrum(i,j,k) = amrex::GpuComplex<amrex::Real>(0.0, 0.0);
                }
            });
        }

        RemapFourierModes(padded_fft_onegrid, original_fft_onegrid,
                          padded_size, original_size);
        original_fft_dist.ParallelCopy(original_fft_onegrid, 0, 0, 1);
        for (amrex::MFIter mfi(original_fft_dist); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.fabbox();
            auto const& spectrum = original_fft_dist.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                spectrum(i,j,k) *= original_scale;
            });
        }
        original_fft.backward(original_fft_dist, filtered_onegrid);
        amrex::MultiFab::Copy(vv_filter, filtered_onegrid, 0, comp, 1, 0);
    }
}


void ProcessDeltaEtaSpectrum(int step, 
                             const amrex::MultiFab& output, 
                             const amrex::Geometry& geom, 
                             amrex::Real kmin, 
                             amrex::Real kmax)
{
    BL_PROFILE("ProcessDeltaEtaSpectrum()");

#if (AMREX_SPACEDIM == 2)
    int const num_comps = 3;
    
    // Extract the real-space delta_eta components (12, 13, 14) from output
    amrex::MultiFab delta_eta_real(output.boxArray(), output.DistributionMap(), num_comps, 0);
    amrex::MultiFab::Copy(delta_eta_real, output, 12, 0, num_comps, 0);

    // Define the metadata needed for IncfloStructFact
    amrex::Vector<std::string> var_names = {"delta_eta_11", "delta_eta_22", "delta_eta_12"};
    amrex::Vector<amrex::Real> var_scaling = {1.0, 1.0, 1.0};
    
    // Pair each variable with itself to compute their auto-spectra
    amrex::Vector<int> pair_a = {0, 1, 2}; 
    amrex::Vector<int> pair_b = {0, 1, 2}; 
    
    // Initialize a local IncfloStructFact object
    IncfloStructFact delta_eta_struct_fact;
    delta_eta_struct_fact.define(output.boxArray(), output.DistributionMap(), geom, 
                                 var_names, var_scaling, pair_a, pair_b, 0);
    
    // Perform the Forward FFT
    delta_eta_struct_fact.sample(delta_eta_real, 1);
    
    // Shift the k=0 mode to the center and compute the magnitude
    // Passing 1 zeroes out the mean (k=0) mode
    delta_eta_struct_fact.callFinalize(1);

    // Format the base file name to include kmin and kmax
    std::ostringstream os;
    os << "Delta_eta_spectrum_" << kmin << "_" << kmax << "_";
    
    // Integrate over shells and write the text file using the custom prefix
    delta_eta_struct_fact.integrateTensorShells(step, os.str());
#else
    amrex::ignore_unused(step, output, geom);
#endif
}


void SpectralWritePlotFile(int step,
                           amrex::Real kmin,
                           amrex::Real kmax,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& velocity,
                           const amrex::MultiFab& velocity_filter,
                           const amrex::MultiFab& vv_filter)
{
    BL_PROFILE_VAR("SpectralWritePlotFile()", SpectralWritePlotFile);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity.nComp() >= AMREX_SPACEDIM,
        "SpectralWritePlotFile: input velocity must have at least AMREX_SPACEDIM components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity_filter.nComp() == 3,
        "SpectralWritePlotFile: velocity_filter must have exactly 3 components");

#if (AMREX_SPACEDIM == 2)
    int constexpr nplot = 15;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vv_filter.nComp() >= 3, "vv_filter must have at least 3 comps in 2D");
#else
    int constexpr nplot = 20;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vv_filter.nComp() >= 6, "vv_filter must have at least 6 comps in 3D");
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

    // amrex::Real sum_tau_S = 0.0;
    // amrex::Real sum_S_sq = 0.0;
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real, amrex::Real> reduce_data(reduce_op);

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

        amrex::Array4<const amrex::Real> const& vv_filt = vv_filter.const_array(mfi);

        // amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        // {
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real, amrex::Real>
        {
            out(i,j,k,0) = vel(i,j,k,0);
            out(i,j,k,1) = vel(i,j,k,1);
#if (AMREX_SPACEDIM == 2)
            out(i,j,k,2) = filt(i,j,k,0);
            out(i,j,k,3) = filt(i,j,k,1);

            amrex::Real const ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
            amrex::Real const vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
            amrex::Real const ux_f = amrex::Real(0.5) * (filt(i+1,j,k,0) - filt(i-1,j,k,0)) * idx;
            amrex::Real const vy_f = amrex::Real(0.5) * (filt(i,j+1,k,1) - filt(i,j-1,k,1)) * idy;

            amrex::Real const vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
            amrex::Real const uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            amrex::Real const vx_f = amrex::Real(0.5) * (filt(i+1,j,k,1) - filt(i-1,j,k,1)) * idx;
            amrex::Real const uy_f = amrex::Real(0.5) * (filt(i,j+1,k,0) - filt(i,j-1,k,0)) * idy;
            out(i,j,k,4) = vx - uy;
            out(i,j,k,5) = vx_f - uy_f;

            out(i,j,k,6) = vv_filt(i,j,k,0); // uu
            out(i,j,k,7) = vv_filt(i,j,k,1); // vv
            out(i,j,k,8) = vv_filt(i,j,k,2); // uv

            // amrex::Real const S_11 = ux;
            amrex::Real const S_11 = ux_f;
            // amrex::Real const S_22 = vy;
            amrex::Real const S_22 = vy_f;
            // amrex::Real const S_12 = amrex::Real(0.5) * (uy + vx);
            amrex::Real const S_12 = amrex::Real(0.5) * (uy_f + vx_f);
            out(i,j,k,9)  = S_11;                               // S_11
            out(i,j,k,10) = S_22;                               // S_22
            
            out(i,j,k,11) =  S_12;    // S_12

            amrex::Real const tau_11 = out(i,j,k,6) - out(i,j,k,2) * out(i,j,k,2);
            amrex::Real const tau_22 = out(i,j,k,7) - out(i,j,k,3) * out(i,j,k,3);
            amrex::Real const tau_12 = out(i,j,k,8) - out(i,j,k,2) * out(i,j,k,3);

            // Calculate half of the trace for the 2D isotropic stress
            amrex::Real const trace_half = amrex::Real(0.5) * (tau_11 + tau_22);

            // Subtract the isotropic part to get the deviatoric stresses
            amrex::Real const tau_11_dev = tau_11 - trace_half;
            amrex::Real const tau_22_dev = tau_22 - trace_half;

            // Deviatoric tau; delta_nu*S is added in below once delta_nu is known
            out(i,j,k,12) = tau_11_dev;
            out(i,j,k,13) = tau_22_dev;
            out(i,j,k,14) = tau_12;
            
            // // Let's compute \tau * S
            // amrex::Real const term1 = (vv_filt(i,j,k,0) - filt(i,j,k,0) * filt(i,j,k,0)) * ux + 
            //                           (vv_filt(i,j,k,1) - filt(i,j,k,1) * filt(i,j,k,1)) * vy + 
            //                           amrex::Real(2.0) *(vv_filt(i,j,k,2) - filt(i,j,k,0) * filt(i,j,k,1)) * S_12;
            
            // Let's compute \tau * S with the deviatoric tau
            amrex::Real const term1 = tau_11_dev * S_11 + 
                                      tau_22_dev * S_22 + 
                                      amrex::Real(2.0) * tau_12 * S_12;
            
            // S^2 as well:
            amrex::Real const term2 = S_11 * S_11 + S_22 * S_22 + amrex::Real(2.0) * S_12 * S_12;
            return {term1, term2};
#else
            out(i,j,k,2) = vel(i,j,k,2);
            out(i,j,k,3) = filt(i,j,k,0);
            out(i,j,k,4) = filt(i,j,k,1);
            out(i,j,k,5) = filt(i,j,k,2);

            amrex::Real const ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
            amrex::Real const vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
            amrex::Real const wz = amrex::Real(0.5) * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

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
            
            out(i,j,k,8)  = vv_filt(i,j,k,0); // uu
            out(i,j,k,9)  = vv_filt(i,j,k,1); // vv
            out(i,j,k,10) = vv_filt(i,j,k,2); // ww
            out(i,j,k,11) = vv_filt(i,j,k,3); // uv
            out(i,j,k,12) = vv_filt(i,j,k,4); // uw
            out(i,j,k,13) = vv_filt(i,j,k,5); // vw

            out(i,j,k,14) = ux;                           // S_11
            out(i,j,k,15) = vy;                           // S_22
            out(i,j,k,16) = wz;                           // S_33
            out(i,j,k,17) = amrex::Real(0.5) * (uy + vx); // S_12
            out(i,j,k,18) = amrex::Real(0.5) * (uz + wx); // S_13
            out(i,j,k,19) = amrex::Real(0.5) * (vz + wy); // S_23
            return {0.0, 0.0};
#endif
        });
    }

#if (AMREX_SPACEDIM == 2)
    amrex::GpuTuple<amrex::Real, amrex::Real> hv = reduce_data.value();
    amrex::Real sum_tau_S = amrex::get<0>(hv);
    amrex::Real sum_S_sq = amrex::get<1>(hv);

    amrex::ParallelDescriptor::ReduceRealSum(sum_tau_S);
    amrex::ParallelDescriptor::ReduceRealSum(sum_S_sq);

    amrex::Real delta_nu = 0.0;
    if (sum_S_sq > 1.0e-20) {
        delta_nu = -(1.0 / (2.0 * sum_S_sq) ) * sum_tau_S;
    }

    // output.setVal(delta_nu, nplot - 1, 1);

    // for (amrex::MFIter mfi(output, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    //     amrex::Box const& bx = mfi.tilebox();
    //     amrex::Array4<amrex::Real> const& out = output.array(mfi);

    //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    //     {
    //         // Reconstruct tau_ij = vv_filt_ij - filt_i * filt_j 
    //         // from the previously saved output components
    //         amrex::Real const tau_11 = out(i,j,k,6) - out(i,j,k,2) * out(i,j,k,2);
    //         amrex::Real const tau_22 = out(i,j,k,7) - out(i,j,k,3) * out(i,j,k,3);
    //         amrex::Real const tau_12 = out(i,j,k,8) - out(i,j,k,2) * out(i,j,k,3);

    //         // Calculate half of the trace for the 2D isotropic stress
    //         amrex::Real const trace_half = amrex::Real(0.5) * (tau_11 + tau_22);

    //         // Subtract the isotropic part to get the deviatoric stresses
    //         amrex::Real const tau_11_dev = tau_11 - trace_half;
    //         amrex::Real const tau_22_dev = tau_22 - trace_half;

    //         // Read the previously saved strain rate tensor components
    //         amrex::Real const S_11 = out(i,j,k,9);
    //         amrex::Real const S_22 = out(i,j,k,10);
    //         amrex::Real const S_12 = out(i,j,k,11);

    //         // Calculate and store delta_eta_ij
    //         // out(i,j,k,12) = tau_11 + amrex::Real(2.0) * delta_nu * S_11;
    //         // out(i,j,k,13) = tau_22 + amrex::Real(2.0) * delta_nu * S_22;
    //         // out(i,j,k,14) = tau_12 + amrex::Real(2.0) * delta_nu * S_12;

    //         // Calculate and store delta_eta_ij using the deviatoric stress
    //         out(i,j,k,12) = tau_11_dev + amrex::Real(2.0) * delta_nu * S_11;
    //         out(i,j,k,13) = tau_22_dev + amrex::Real(2.0) * delta_nu * S_22;
    //         out(i,j,k,14) = tau_12 + amrex::Real(2.0) * delta_nu * S_12;
    //     });
    // }
    // delta_eta_ij = tau_ij_dev + 2 * delta_nu * S_ij   (comps 12..14 += 2*dnu * comps 9..11)
    amrex::MultiFab::Saxpy(output, amrex::Real(2.0) * delta_nu, output, 9, 12, 3, 0);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string const dNufilename = Delta_nu_FileName(step);

        static std::set<std::string> opened;   // truncate once per run, then append
        bool const first = opened.insert(dNufilename).second;
        std::ofstream ofs(dNufilename, first ? std::ios::trunc : std::ios::app);
        if (!ofs) {
            amrex::Abort("Could not open " + dNufilename);
        }
        if (first) { ofs << "# step  k_min  k_cutoff  delta_nu  sum_tau_S  sum_S_sq\n"; }
        ofs.precision(17);
        ofs << step << " " << kmin << " " << kmax << " " << delta_nu << " "
            << sum_tau_S << " " << sum_S_sq << "\n";
    }
            
    amrex::Vector<std::string> varNames{
        "velx", "vely",
        "velx_filter", "vely_filter",
        "vort", "vort_filter",
        "uu_filter", "vv_filter", "uv_filter",
        "S11", "S22", "S12", 
        "delta_eta_11", "delta_eta_22", "delta_eta_12"};

    ProcessDeltaEtaSpectrum(step, output, geom, kmin, kmax);
#else
    amrex::Vector<std::string> varNames{
        "velx", "vely", "velz",
        "velx_filter", "vely_filter", "velz_filter",
        "vort", "vort_filter",
        "uu_filter", "vv_filter", "ww_filter",
        "uv_filter", "uw_filter", "vw_filter",
        "S11", "S22", "S33", "S12", "S13", "S23"};
#endif

    std::string const plotfilename = FilterPlotFileName(step, kmin, kmax);
    amrex::Print() << "Writing filtered velocity plotfile " << plotfilename << "\n";
    amrex::WriteSingleLevelPlotfile(plotfilename, output, varNames, geom, 0.0, step);

    
}


void SpectralWriteFourierPlotFile(int step,
                                  amrex::Real kmin,
                                  amrex::Real kmax,
                                  const amrex::Geometry& geom,
                                  const amrex::MultiFab& velocity_filter,
                                  const amrex::MultiFab& vv_filter)
{
    BL_PROFILE_VAR("SpectralWriteFourierPlotFile()", SpectralWriteFourierPlotFile);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        velocity_filter.nComp() == 3,
        "SpectralWriteFourierPlotFile: velocity_filter must have exactly 3 components");

#if (AMREX_SPACEDIM == 2)
    int constexpr num_vv_comps = 3; 
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vv_filter.nComp() >= 3, "vv_filter must have at least 3 comps in 2D");
#else
    int constexpr num_vv_comps = 6;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vv_filter.nComp() >= 6, "vv_filter must have at least 6 comps in 3D");
#endif

    // int constexpr nfields = AMREX_SPACEDIM + 1;
    int constexpr nfields = AMREX_SPACEDIM + 1 + num_vv_comps;

    amrex::MultiFab velocity_filter_g(velocity_filter.boxArray(), velocity_filter.DistributionMap(), 3, 1);
    FillVelocityWithGhosts(velocity_filter, velocity_filter_g, geom);

    amrex::MultiFab fields(velocity_filter.boxArray(), velocity_filter.DistributionMap(), nfields, 0);
    fields.setVal(0.0);
    amrex::MultiFab::Copy(fields, velocity_filter, 0, 0, AMREX_SPACEDIM, 0);

    amrex::MultiFab vort_filter(fields, amrex::make_alias, AMREX_SPACEDIM, 1);
    ComputeVorticityField(velocity_filter_g, vort_filter, geom);

    amrex::MultiFab::Copy(fields, vv_filter, 0, AMREX_SPACEDIM + 1, num_vv_comps, 0);

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
        "velx_filter", "vely_filter", "vort_filter",
        "uu_filter", "vv_filter", "uv_filter"};
#else
    amrex::Vector<std::string> baseNames{
        "velx_filter", "vely_filter", "velz_filter", "vort_filter",
        "uu_filter", "vv_filter", "ww_filter",
        "uv_filter", "uw_filter", "vw_filter"};
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
