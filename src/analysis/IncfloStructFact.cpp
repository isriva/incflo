#include "IncfloStructFact.H"

#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

using namespace amrex;

namespace {

void
goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void
sqrt_mf(MultiFab& mf)
{
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);
        int const ncomp = mf.nComp();

        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            arr(i,j,k,n) = std::sqrt(arr(i,j,k,n));
        });
    }
}

}

void
IncfloStructFact::define(BoxArray const& ba,
                         DistributionMapping const& dm,
                         Geometry const& geom,
                         Vector<std::string> const& var_names,
                         Vector<Real> const& var_scaling,
                         int verbosity)
{
    int const nvar = static_cast<int>(var_names.size());
    Vector<int> pair_a(nvar * (nvar + 1) / 2);
    Vector<int> pair_b(nvar * (nvar + 1) / 2);

    int counter = 0;
    for (int i = 0; i < nvar; ++i) {
        for (int j = i; j < nvar; ++j) {
            pair_a[counter] = j;
            pair_b[counter] = i;
            ++counter;
        }
    }

    define(ba, dm, geom, var_names, var_scaling, pair_a, pair_b, verbosity);
}

void
IncfloStructFact::define(BoxArray const& ba,
                         DistributionMapping const& dm,
                         Geometry const& geom,
                         Vector<std::string> const& var_names,
                         Vector<Real> const& var_scaling,
                         Vector<int> const& pair_a,
                         Vector<int> const& pair_b,
                         int verbosity)
{
    BL_PROFILE("IncfloStructFact::define()");

    m_verbosity = verbosity;

    if (pair_a.size() != pair_b.size()) {
        Abort("IncfloStructFact::define: pair vectors must have the same length");
    }

    m_nvar = static_cast<int>(var_names.size());
    m_ncov = static_cast<int>(pair_a.size());

    if (m_ncov != static_cast<int>(var_scaling.size())) {
        Abort("IncfloStructFact::define: covariance count does not match scaling count");
    }

    m_scaling.resize(m_ncov);
    m_pair_a.resize(m_ncov);
    m_pair_b.resize(m_ncov);

    for (int n = 0; n < m_ncov; ++n) {
        if (var_scaling[n] == Real(0.0)) {
            Abort("IncfloStructFact::define: zero covariance scaling");
        }
        m_scaling[n] = Real(1.0) / var_scaling[n];
        m_pair_a[n] = pair_a[n];
        m_pair_b[n] = pair_b[n];
        if (m_pair_a[n] < 0 || m_pair_a[n] >= m_nvar ||
            m_pair_b[n] < 0 || m_pair_b[n] >= m_nvar) {
            Abort("IncfloStructFact::define: covariance pair index is out of range");
        }
    }

    Vector<int> var_unique_temp;
    var_unique_temp.reserve(2 * m_ncov);
    for (int n = 0; n < m_ncov; ++n) {
        var_unique_temp.push_back(m_pair_a[n]);
        var_unique_temp.push_back(m_pair_b[n]);
    }
    std::sort(var_unique_temp.begin(), var_unique_temp.end());
    var_unique_temp.erase(std::unique(var_unique_temp.begin(), var_unique_temp.end()),
                          var_unique_temp.end());

    m_nvar_unique = static_cast<int>(var_unique_temp.size());
    m_var_unique = var_unique_temp;

    if (m_verbosity > 0) {
        for (int n = 0; n < m_ncov; ++n) {
            Print() << "IncfloStructFact pair " << n << " = "
                    << m_pair_b[n] << ", " << m_pair_a[n] << "\n";
        }
        Print() << "IncfloStructFact numPairs = " << m_ncov << "\n";
    }

    m_cov_real.define(ba, dm, m_ncov, 0);
    m_cov_imag.define(ba, dm, m_ncov, 0);
    m_cov_mag.define(ba, dm, m_ncov, 0);
    reset();

    m_cov_names.resize(m_ncov);
    for (int n = 0; n < m_ncov; ++n) {
        m_cov_names[n] = "struct_fact_" + var_names[m_pair_b[n]] + "_" +
                         var_names[m_pair_a[n]];
    }

    Box const domain = geom.Domain();
    m_ncell = domain.size();

    Vector<Real> kspace_lo(AMREX_SPACEDIM);
    Vector<Real> kspace_hi(AMREX_SPACEDIM);

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (domain.length(d) % 2 == 0) {
            kspace_lo[d] = -domain.length(d) / Real(2.0) - Real(0.5);
            kspace_hi[d] =  domain.length(d) / Real(2.0) - Real(0.5);
        } else {
            kspace_lo[d] = -domain.length(d) / Real(2.0);
            kspace_hi[d] =  domain.length(d) / Real(2.0);
        }
    }

    RealBox kspace({AMREX_D_DECL(kspace_lo[0], kspace_lo[1], kspace_lo[2])},
                   {AMREX_D_DECL(kspace_hi[0], kspace_hi[1], kspace_hi[2])});

    Vector<int> is_periodic(AMREX_SPACEDIM, 1);
    m_kspace_geom.define(domain, &kspace, CoordSys::cartesian, is_periodic.data());
}

void
IncfloStructFact::sample(MultiFab const& variables, int reset_sample)
{
    BL_PROFILE("IncfloStructFact::sample()");

    BoxArray const& ba = variables.boxArray();
    DistributionMapping const& dm = variables.DistributionMap();

    MultiFab variables_dft_real(ba, dm, m_nvar, 0);
    MultiFab variables_dft_imag(ba, dm, m_nvar, 0);

    computeFFT(variables, variables_dft_real, variables_dft_imag);

    MultiFab cov_temp(ba, dm, 1, 0);
    MultiFab cov_temp2(m_cov_real.boxArray(), m_cov_real.DistributionMap(), 1, 0);

    for (int n = 0; n < m_ncov; ++n) {
        int const i = m_pair_a[n];
        int const j = m_pair_b[n];

        cov_temp.setVal(0.0);
        MultiFab::AddProduct(cov_temp, variables_dft_real, i, variables_dft_real, j, 0, 1, 0);
        MultiFab::AddProduct(cov_temp, variables_dft_imag, i, variables_dft_imag, j, 0, 1, 0);
        cov_temp2.ParallelCopy(cov_temp, 0, 0, 1);

        if (reset_sample == 1) {
            MultiFab::Copy(m_cov_real, cov_temp2, 0, n, 1, 0);
        } else {
            MultiFab::Add(m_cov_real, cov_temp2, 0, n, 1, 0);
        }

        cov_temp.setVal(0.0);
        MultiFab::AddProduct(cov_temp, variables_dft_imag, i, variables_dft_real, j, 0, 1, 0);
        cov_temp.mult(-1.0, 0);
        MultiFab::AddProduct(cov_temp, variables_dft_real, i, variables_dft_imag, j, 0, 1, 0);
        cov_temp2.ParallelCopy(cov_temp, 0, 0, 1);

        if (reset_sample == 1) {
            MultiFab::Copy(m_cov_imag, cov_temp2, 0, n, 1, 0);
        } else {
            MultiFab::Add(m_cov_imag, cov_temp2, 0, n, 1, 0);
        }
    }

    m_nsamples = (reset_sample == 1) ? 1 : m_nsamples + 1;
}

void
IncfloStructFact::reset()
{
    BL_PROFILE("IncfloStructFact::reset()");

    m_cov_real.setVal(0.0);
    m_cov_imag.setVal(0.0);
    m_cov_mag.setVal(0.0);
    m_nsamples = 0;
}

void
IncfloStructFact::computeFFT(MultiFab const& variables,
                             MultiFab& variables_dft_real,
                             MultiFab& variables_dft_imag,
                             bool unpack) const
{
    BL_PROFILE("IncfloStructFact::computeFFT()");

    Box const domain = variables.boxArray().minimalBox();
    bool const chopped_in_x = domain.length(0) > 1;
    bool const chopped_in_y = !chopped_in_x && domain.length(1) > 1;
#if (AMREX_SPACEDIM == 3)
    bool const chopped_in_z = !chopped_in_x && !chopped_in_y && domain.length(2) > 1;
#else
    bool const chopped_in_z = false;
#endif

    if (!chopped_in_x && !chopped_in_y && !chopped_in_z) {
        Abort("IncfloStructFact::computeFFT: cannot transform a one-cell domain");
    }

    Long npts = AMREX_D_TERM(Long(domain.length(0)),
                             * Long(domain.length(1)),
                             * Long(domain.length(2)));
    Real const sqrtnpts = std::sqrt(Real(npts));

    BoxArray const ba = variables.boxArray();
    DistributionMapping const dm = variables.DistributionMap();

    MultiFab phi(ba, dm, 1, 0);

    BoxArray ba_onegrid(domain);
    DistributionMapping dm_onegrid(ba_onegrid);

    FFT::R2C my_fft(domain);
    auto const& [ba_fft, dm_fft] = my_fft.getSpectralDataLayout();
    FabArray<BaseFab<GpuComplex<Real>>> phi_fft(ba_fft, dm_fft, 1, 0);

    Box const domain_fft = ba_fft.minimalBox();
    BoxArray ba_fft_onegrid(domain_fft);
    FabArray<BaseFab<GpuComplex<Real>>> phi_fft_onegrid(ba_fft_onegrid, dm_onegrid, 1, 0);

    MultiFab variables_dft_real_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    MultiFab variables_dft_imag_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    for (int comp = 0; comp < m_nvar; ++comp) {
        bool comp_fft = false;
        for (int i = 0; i < m_nvar_unique; ++i) {
            if (comp == m_var_unique[i]) {
                comp_fft = true;
                break;
            }
        }
        if (!comp_fft) {
            continue;
        }

        MultiFab::Copy(phi, variables, comp, 0, 1, 0);
        my_fft.forward(phi, phi_fft);
        phi_fft_onegrid.ParallelCopy(phi_fft, 0, 0, 1);

        for (MFIter mfi(variables_dft_real_onegrid); mfi.isValid(); ++mfi) {
            Box const bx = mfi.fabbox();
            Array4<GpuComplex<Real>> const spectral = phi_fft_onegrid.array(mfi);
            Array4<Real> const realpart = variables_dft_real_onegrid.array(mfi);
            Array4<Real> const imagpart = variables_dft_imag_onegrid.array(mfi);

            if (chopped_in_x) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                        realpart(i,j,k) = unpack ?  spectral(iloc,jloc,kloc).real() : Real(0.0);
                        imagpart(i,j,k) = unpack ? -spectral(iloc,jloc,kloc).imag() : Real(0.0);
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
            }

            if (chopped_in_y) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                        realpart(i,j,k) = unpack ?  spectral(iloc,jloc,kloc).real() : Real(0.0);
                        imagpart(i,j,k) = unpack ? -spectral(iloc,jloc,kloc).imag() : Real(0.0);
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
            }

#if (AMREX_SPACEDIM == 3)
            if (chopped_in_z) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k <= bx.length(2) / 2) {
                        realpart(i,j,k) = spectral(i,j,k).real();
                        imagpart(i,j,k) = spectral(i,j,k).imag();
                    } else {
                        int const iloc = (i == 0) ? 0 : bx.length(0) - i;
                        int const jloc = (j == 0) ? 0 : bx.length(1) - j;
                        int const kloc = bx.length(2) - k;
                        realpart(i,j,k) = unpack ?  spectral(iloc,jloc,kloc).real() : Real(0.0);
                        imagpart(i,j,k) = unpack ? -spectral(iloc,jloc,kloc).imag() : Real(0.0);
                    }
                    realpart(i,j,k) /= sqrtnpts;
                    imagpart(i,j,k) /= sqrtnpts;
                });
            }
#endif
        }

        variables_dft_real.ParallelCopy(variables_dft_real_onegrid, 0, comp, 1);
        variables_dft_imag.ParallelCopy(variables_dft_imag_onegrid, 0, comp, 1);
    }
}

void
IncfloStructFact::writePlotFile(int step, Real time,
                                std::string const& plotfile_base,
                                int zero_mode)
{
    BL_PROFILE("IncfloStructFact::writePlotFile()");

    if (m_nsamples <= 0) {
        return;
    }

    BoxArray const& ba = m_cov_mag.boxArray();
    DistributionMapping const& dm = m_cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, m_ncov, 0);
    MultiFab cov_imag_temp(ba, dm, m_ncov, 0);
    MultiFab::Copy(cov_real_temp, m_cov_real, 0, 0, m_ncov, 0);
    MultiFab::Copy(cov_imag_temp, m_cov_imag, 0, 0, m_ncov, 0);

    finalize(cov_real_temp, cov_imag_temp, zero_mode);

    {
        std::string const plotfilename = Concatenate(plotfile_base + "_mag", step, 9);
        MultiFab plotfile(ba, dm, m_ncov, 0);
        MultiFab::Copy(plotfile, m_cov_mag, 0, 0, m_ncov, 0);
        WriteSingleLevelPlotfile(plotfilename, plotfile, m_cov_names, m_kspace_geom, time, step);
    }

    {
        std::string const plotfilename = Concatenate(plotfile_base + "_real_imag", step, 9);
        MultiFab plotfile(ba, dm, 2 * m_ncov, 0);
        Vector<std::string> var_names(2 * m_ncov);

        for (int n = 0; n < m_ncov; ++n) {
            var_names[n] = m_cov_names[n] + "_real";
            var_names[n + m_ncov] = m_cov_names[n] + "_imag";
        }

        MultiFab::Copy(plotfile, cov_real_temp, 0, 0, m_ncov, 0);
        MultiFab::Copy(plotfile, cov_imag_temp, 0, m_ncov, m_ncov, 0);
        WriteSingleLevelPlotfile(plotfilename, plotfile, var_names, m_kspace_geom, time, step);
    }
}

void
IncfloStructFact::finalize(MultiFab& cov_real_in, MultiFab& cov_imag_in, int zero_mode)
{
    BL_PROFILE("IncfloStructFact::finalize()");

    if (m_nsamples <= 0) {
        Abort("IncfloStructFact::finalize: no samples available");
    }

    Real const nsamples_inv = Real(1.0) / Real(m_nsamples);

    shiftFFT(cov_real_in, zero_mode);
    shiftFFT(cov_imag_in, zero_mode);

    cov_real_in.mult(nsamples_inv);
    cov_imag_in.mult(nsamples_inv);
    for (int n = 0; n < m_ncov; ++n) {
        cov_real_in.mult(m_scaling[n], n, 1);
        cov_imag_in.mult(m_scaling[n], n, 1);
    }

    m_cov_mag.setVal(0.0);
    MultiFab::AddProduct(m_cov_mag, cov_real_in, 0, cov_real_in, 0, 0, m_ncov, 0);
    MultiFab::AddProduct(m_cov_mag, cov_imag_in, 0, cov_imag_in, 0, 0, m_ncov, 0);
    sqrt_mf(m_cov_mag);
}

void
IncfloStructFact::callFinalize(int zero_mode)
{
    BL_PROFILE("IncfloStructFact::callFinalize()");

    BoxArray const& ba = m_cov_mag.boxArray();
    DistributionMapping const& dm = m_cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, m_ncov, 0);
    MultiFab cov_imag_temp(ba, dm, m_ncov, 0);
    MultiFab::Copy(cov_real_temp, m_cov_real, 0, 0, m_ncov, 0);
    MultiFab::Copy(cov_imag_temp, m_cov_imag, 0, 0, m_ncov, 0);

    finalize(cov_real_temp, cov_imag_temp, zero_mode);
}

void
IncfloStructFact::shiftFFT(MultiFab& dft_out, int zero_mode) const
{
    BL_PROFILE("IncfloStructFact::shiftFFT()");

    Box const domain = dft_out.boxArray().minimalBox();
    BoxArray ba_onegrid(domain);
    DistributionMapping dm_onegrid(ba_onegrid);

    MultiFab dft_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    MultiFab dft_onegrid_temp(ba_onegrid, dm_onegrid, 1, 0);

    for (int n = 0; n < m_ncov; ++n) {
        dft_onegrid_temp.ParallelCopy(dft_out, n, 0, 1);

        if (zero_mode == 1) {
            for (MFIter mfi(dft_onegrid_temp); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const dft_temp = dft_onegrid_temp.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i == 0 && j == 0 && k == 0) {
                        dft_temp(i,j,k) = Real(0.0);
                    }
                });
            }
        }

        for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            Array4<Real> const dft = dft_onegrid.array(mfi);
            Array4<Real const> const dft_temp = dft_onegrid_temp.const_array(mfi);

            int const nx = bx.length(0);
            int const nxh = nx / 2;
            int const ny = bx.length(1);
            int const nyh = ny / 2;
#if (AMREX_SPACEDIM == 3)
            int const nz = bx.length(2);
            int const nzh = nz / 2;
#endif

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int const ip = (i + nxh) % nx;
                int const jp = (j + nyh) % ny;
                int kp = 0;
#if (AMREX_SPACEDIM == 3)
                kp = (k + nzh) % nz;
#endif
                dft(ip,jp,kp) = dft_temp(i,j,k);
            });
        }

        dft_out.ParallelCopy(dft_onegrid, 0, n, 1);
    }
}

void
IncfloStructFact::integrateShells(int step, std::string const& output_base) const
{
    BL_PROFILE("IncfloStructFact::integrateShells()");

    if (m_nsamples <= 0) {
        return;
    }

    GpuArray<int, AMREX_SPACEDIM> center;
    int npts = std::numeric_limits<int>::max();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        center[d] = m_ncell[d] / 2;
        npts = std::min(npts, center[d]);
    }

    if (npts <= 1) {
        return;
    }

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int> phicnt_device(npts);
    Gpu::HostVector<Real> phisum_host(npts);

    Real* phisum_ptr = phisum_device.dataPtr();
    int* phicnt_ptr = phicnt_device.dataPtr();

    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        phisum_ptr[d] = Real(0.0);
        phicnt_ptr[d] = 0;
    });

    int const ncomp_sum = std::min(AMREX_SPACEDIM, m_ncov);
    for (MFIter mfi(m_cov_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real const> const cov = m_cov_mag.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int const ilen = amrex::Math::abs(i - center[0]);
            int const jlen = amrex::Math::abs(j - center[1]);
            int const klen = (AMREX_SPACEDIM == 3) ? amrex::Math::abs(k - center[2]) : 0;

            Real dist = std::sqrt(Real(ilen*ilen + jlen*jlen + klen*klen));
            if (dist <= Real(npts) - Real(0.5)) {
                int const shell = int(dist + Real(0.5));
                for (int n = 0; n < ncomp_sum; ++n) {
                    HostDevice::Atomic::Add(&(phisum_ptr[shell]), cov(i,j,k,n));
                }
                HostDevice::Atomic::Add(&(phicnt_ptr[shell]), 1);
            }
        });
    }

    for (int d = 1; d < npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }

    Real const pi = Real(4.0) * std::atan(Real(1.0));
    Real const dk = Real(1.0);

#if (AMREX_SPACEDIM == 2)
    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0 && phicnt_ptr[d] > 0) {
            phisum_ptr[d] *= Real(2.0) * pi * (d * dk + Real(0.5) * dk * dk) /
                             phicnt_ptr[d];
        }
    });
#else
    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0 && phicnt_ptr[d] > 0) {
            phisum_ptr[d] *= Real(4.0) * pi *
                             (d * d * dk + dk * dk * dk / Real(12.0)) /
                             phicnt_ptr[d];
        }
    });
#endif

    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(),
              phisum_host.begin());

    if (ParallelDescriptor::IOProcessor()) {
        std::string filename = Concatenate(output_base, step, 7) + ".txt";
        std::ofstream turb(filename);
        turb.precision(17);
        for (int d = 1; d < npts; ++d) {
            turb << d << " " << phisum_host[d] << "\n";
        }
    }
}


void
IncfloStructFact::integrateTensorShells(int step, std::string const& output_base) const
{
    BL_PROFILE("IncfloStructFact::integrateTensorShells()");

    if (m_nsamples <= 0) {
        return;
    }

    GpuArray<int, AMREX_SPACEDIM> center;
    int npts = std::numeric_limits<int>::max();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        center[d] = m_ncell[d] / 2;
        npts = std::min(npts, center[d]);
    }

    if (npts <= 1) {
        return;
    }

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int> phicnt_device(npts);
    Gpu::HostVector<Real> phisum_host(npts);

    Real* phisum_ptr = phisum_device.dataPtr();
    int* phicnt_ptr = phicnt_device.dataPtr();

    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        phisum_ptr[d] = Real(0.0);
        phicnt_ptr[d] = 0;
    });

    int const ncomp_sum = m_ncov;
    for (MFIter mfi(m_cov_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real const> const cov = m_cov_mag.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int const ilen = amrex::Math::abs(i - center[0]);
            int const jlen = amrex::Math::abs(j - center[1]);
            int const klen = (AMREX_SPACEDIM == 3) ? amrex::Math::abs(k - center[2]) : 0;

            Real dist = std::sqrt(Real(ilen*ilen + jlen*jlen + klen*klen));
            if (dist <= Real(npts) - Real(0.5)) {
                int const shell = int(dist + Real(0.5));
                // for (int n = 0; n < ncomp_sum; ++n) {
                //     HostDevice::Atomic::Add(&(phisum_ptr[shell]), cov(i,j,k,n));
                // }
                // HostDevice::Atomic::Add(&(phicnt_ptr[shell]), 1);
                Real tensor_energy = Real(0.0);
                for (int n = 0; n < ncomp_sum; ++n) {
                    Real amp = cov(i,j,k,n);
                    Real weight = (n == 2) ? Real(2.0) : Real(1.0);
                    tensor_energy += weight * (amp * amp);
                }
                
                HostDevice::Atomic::Add(&(phisum_ptr[shell]), tensor_energy);
                HostDevice::Atomic::Add(&(phicnt_ptr[shell]), 1);
            }
        });
    }

    for (int d = 1; d < npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }

    Real const pi = Real(4.0) * std::atan(Real(1.0));
    Real const dk = Real(1.0);

#if (AMREX_SPACEDIM == 2)
    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0 && phicnt_ptr[d] > 0) {
            phisum_ptr[d] *= Real(2.0) * pi * (d * dk + Real(0.5) * dk * dk) /
                             phicnt_ptr[d];
        }
    });
#else
    ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0 && phicnt_ptr[d] > 0) {
            phisum_ptr[d] *= Real(4.0) * pi *
                             (d * d * dk + dk * dk * dk / Real(12.0)) /
                             phicnt_ptr[d];
        }
    });
#endif

    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(),
              phisum_host.begin());

    if (ParallelDescriptor::IOProcessor()) {
        // std::string filename = Concatenate(output_base, step, 7) + ".txt";
        std::string filename = output_base + std::to_string(step) + ".txt";
        std::ofstream turb(filename);
        turb.precision(17);
        for (int d = 1; d < npts; ++d) {
            turb << d << " " << phisum_host[d] << "\n";
        }
    }
}


void
IncfloStructFact::writeCheckpoint(std::string const& checkpoint_dir) const
{
    BL_PROFILE("IncfloStructFact::writeCheckpoint()");

    PreBuildDirectorHierarchy(checkpoint_dir, "Level_", 1, true);

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream header;
        header.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string const header_name = checkpoint_dir + "/Header";
        header.open(header_name.c_str(), std::ofstream::out |
                    std::ofstream::trunc | std::ofstream::binary);
        if (!header.good()) {
            FileOpenFailed(header_name);
        }

        header.precision(17);
        header << "incflo structure factor checkpoint\n";
        header << 2 << "\n";
        header << m_nvar << "\n";
        header << m_nvar_unique << "\n";
        header << m_ncov << "\n";
        header << m_nsamples << "\n";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            header << m_ncell[d] << " ";
        }
        header << "\n";
        for (int n = 0; n < m_ncov; ++n) {
            header << m_cov_names[n] << "\n";
        }
        for (int n = 0; n < m_ncov; ++n) {
            header << m_pair_a[n] << " " << m_pair_b[n] << " "
                   << m_scaling[n] << "\n";
        }
    }

    VisMF::Write(m_cov_real, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_real"));
    VisMF::Write(m_cov_imag, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_imag"));
    VisMF::Write(m_cov_mag, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_mag"));
}

void
IncfloStructFact::readCheckpoint(std::string const& checkpoint_dir)
{
    BL_PROFILE("IncfloStructFact::readCheckpoint()");

    std::string line;
    int version = 0;
    int nvar = 0;
    int nvar_unique = 0;
    int ncov = 0;
    int nsamples = 0;
    IntVect ncell{AMREX_D_DECL(1, 1, 1)};
    Vector<std::string> cov_names;
    Vector<int> pair_a;
    Vector<int> pair_b;
    Vector<Real> scaling;

    {
        std::string const file = checkpoint_dir + "/Header";
        Vector<char> file_char_ptr;
        ParallelDescriptor::ReadAndBcastFile(file, file_char_ptr);
        std::string file_string(file_char_ptr.dataPtr());
        std::istringstream is(file_string, std::istringstream::in);

        std::getline(is, line);
        is >> version; goto_next_line(is);
        is >> nvar; goto_next_line(is);
        is >> nvar_unique; goto_next_line(is);
        is >> ncov; goto_next_line(is);
        is >> nsamples; goto_next_line(is);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            is >> ncell[d];
        }
        goto_next_line(is);

        cov_names.resize(ncov);
        pair_a.resize(ncov);
        pair_b.resize(ncov);
        scaling.resize(ncov);
        for (int n = 0; n < ncov; ++n) {
            is >> cov_names[n];
            goto_next_line(is);
        }
        for (int n = 0; n < ncov; ++n) {
            is >> pair_a[n] >> pair_b[n] >> scaling[n];
            goto_next_line(is);
        }
    }

    if (version != 2 || nvar != m_nvar || nvar_unique != m_nvar_unique ||
        ncov != m_ncov) {
        Abort("IncfloStructFact::readCheckpoint: checkpoint metadata does not match current analysis configuration");
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (ncell[d] != m_ncell[d]) {
            Abort("IncfloStructFact::readCheckpoint: checkpoint domain does not match current analysis configuration");
        }
    }
    for (int n = 0; n < m_ncov; ++n) {
        Real const scale_tol = Real(64.0) * std::numeric_limits<Real>::epsilon() *
            std::max(Real(1.0), std::max(std::abs(scaling[n]), std::abs(m_scaling[n])));
        if (cov_names[n] != m_cov_names[n] ||
            pair_a[n] != m_pair_a[n] || pair_b[n] != m_pair_b[n] ||
            std::abs(scaling[n] - m_scaling[n]) > scale_tol) {
            Abort("IncfloStructFact::readCheckpoint: covariance metadata does not match current analysis configuration");
        }
    }

    m_nsamples = nsamples;

    VisMF::Read(m_cov_real, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_real"));
    VisMF::Read(m_cov_imag, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_imag"));
    VisMF::Read(m_cov_mag, MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_mag"));
}

