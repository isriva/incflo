#include <incflo.H>

#ifdef INCFLO_USE_FFT
#include "analysis/IncfloStructFact.H"
#include <AMReX_VisMF.H>
#endif

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <fstream>
#include <memory>

using namespace amrex;

namespace {

Vector<std::string>
velocity_names()
{
#if (AMREX_SPACEDIM == 3)
    return {"velx", "vely", "velz"};
#else
    return {"velx", "vely"};
#endif
}

Real
velocity_cell_volume(Geometry const& geom, Real cell_depth)
{
    auto const dx = geom.CellSizeArray();
    return AMREX_D_TERM(dx[0], * dx[1], * dx[2]) *
#if (AMREX_SPACEDIM == 2)
           cell_depth;
#else
           Real(1.0);
#endif
}


#ifdef INCFLO_USE_FFT
bool
struct_fact_checkpoint_exists(std::string const& checkpoint_dir)
{
    std::ifstream header(checkpoint_dir + "/Header");
    if (!header.good()) {
        return false;
    }

    std::string const cov_real =
        MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_real");
    std::string const cov_imag =
        MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_imag");
    std::string const cov_mag =
        MultiFabFileFullPrefix(0, checkpoint_dir, "Level_", "cov_mag");
    return VisMF::Exist(cov_real) && VisMF::Exist(cov_imag) &&
           VisMF::Exist(cov_mag);
}
#endif

}

bool
incflo::structFactEnabled () const
{
    return m_struct_fact_int > 0 || m_turb_spectrum_int > 0;
}

bool
incflo::structFactSampleNow () const
{
    if (m_struct_fact_int <= 0) {
        return false;
    }
    if (m_nstep == 0) {
        return m_struct_fact_skip == 0;
    }
    return m_nstep > m_struct_fact_skip &&
           (m_nstep - m_struct_fact_skip) % m_struct_fact_int == 0;
}

bool
incflo::structFactWriteNow () const
{
    if (m_struct_fact_int <= 0 || m_struct_fact_write_int <= 0) {
        return false;
    }
    if (m_nstep == 0) {
        return m_struct_fact_skip == 0;
    }
    return m_nstep > m_struct_fact_skip &&
           m_nstep % m_struct_fact_write_int == 0;
}

bool
incflo::turbSpectrumWriteNow () const
{
    return m_turb_spectrum_int > 0 &&
           m_nstep >= 0 &&
           m_nstep % m_turb_spectrum_int == 0;
}

void
incflo::InitStructFact (bool restarted)
{
    if (!structFactEnabled()) {
        return;
    }

#if !defined(INCFLO_USE_FFT) || !defined(AMREX_USE_FFT)
    amrex::Abort("incflo structure-factor or turbulent-spectrum analysis requires an FFT-enabled incflo build");
#else
    if (finest_level != 0 || max_level > 0) {
        amrex::Abort("incflo structure-factor analysis requires max_level = 0");
    }

    Box const& domain = Geom(0).Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (!Geom(0).isPeriodic(dir)) {
            amrex::Abort("incflo structure-factor analysis requires a fully periodic domain");
        }
        if (domain.smallEnd(dir) != 0) {
            amrex::Abort("incflo structure-factor analysis requires domain smallEnd = 0");
        }
    }

#ifdef AMREX_USE_EB
    if (!EBFactory(0).isAllRegular()) {
        amrex::Abort("incflo structure-factor analysis requires all-regular EB geometry");
    }
#endif

    if (m_struct_fact_cell_depth <= Real(0.0)) {
        amrex::Abort("incflo structure-factor analysis requires struct_fact_cell_depth > 0");
    }

    Vector<std::string> const var_names = velocity_names();
    Real const dvol = velocity_cell_volume(Geom(0), m_struct_fact_cell_depth);

    if (m_struct_fact_int > 0) {
        int const npairs = AMREX_SPACEDIM * (AMREX_SPACEDIM + 1) / 2;
        Vector<Real> scaling(npairs, Real(1.0) / dvol);

        m_struct_fact = std::make_unique<IncfloStructFact>();
        m_struct_fact->define(m_leveldata[0]->velocity.boxArray(),
                              m_leveldata[0]->velocity.DistributionMap(),
                              Geom(0), var_names, scaling,
                              m_struct_fact_verbosity);

        if (restarted) {
            std::string const checkpoint_dir = m_restart_file + "/StructFact";
            if (struct_fact_checkpoint_exists(checkpoint_dir)) {
                m_struct_fact->readCheckpoint(checkpoint_dir);
                m_struct_fact_needs_write = m_struct_fact->samples() > 0;
            } else {
                amrex::Warning("StructFact restart state is missing or incomplete; starting a fresh structure-factor average");
            }
        }
    }

    if (m_turb_spectrum_int > 0) {
        Vector<int> pair_a(AMREX_SPACEDIM);
        Vector<int> pair_b(AMREX_SPACEDIM);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            pair_a[dir] = dir;
            pair_b[dir] = dir;
        }
        Vector<Real> scaling(AMREX_SPACEDIM, Real(1.0) / dvol);

        m_turb_spectrum = std::make_unique<IncfloStructFact>();
        m_turb_spectrum->define(m_leveldata[0]->velocity.boxArray(),
                                m_leveldata[0]->velocity.DistributionMap(),
                                Geom(0), var_names, scaling, pair_a, pair_b,
                                m_struct_fact_verbosity);
    }
#endif
}

void
incflo::SampleStructFact ()
{
#ifdef INCFLO_USE_FFT
    if (m_struct_fact && structFactSampleNow()) {
        m_struct_fact->sample(m_leveldata[0]->velocity);
        m_last_struct_fact_sample = m_nstep;
        m_struct_fact_needs_write = true;
    }
#endif
}

void
incflo::WriteStructFact (bool force)
{
#ifdef INCFLO_USE_FFT
    if (m_struct_fact &&
        m_struct_fact->samples() > 0 &&
        m_struct_fact_needs_write &&
        (force || structFactWriteNow()) &&
        m_last_struct_fact_write != m_nstep) {
        m_struct_fact->writePlotFile(m_nstep, m_cur_time, m_struct_fact_file,
                                     m_struct_fact_zero_mode);
        m_last_struct_fact_write = m_nstep;
        m_struct_fact_needs_write = false;
    }
#else
    amrex::ignore_unused(force);
#endif
}

void
incflo::WriteTurbulentSpectrum ()
{
#ifdef INCFLO_USE_FFT
    if (m_turb_spectrum && turbSpectrumWriteNow()) {
        m_turb_spectrum->sample(m_leveldata[0]->velocity, 1);
        m_turb_spectrum->callFinalize(m_struct_fact_zero_mode);
        m_turb_spectrum->integrateShells(m_nstep, m_turb_spectrum_file);
    }
#endif
}

void
incflo::WriteStructFactCheckpoint (std::string const& checkpointname) const
{
#ifdef INCFLO_USE_FFT
    if (m_struct_fact) {
        m_struct_fact->writeCheckpoint(checkpointname + "/StructFact");
    }
#else
    amrex::ignore_unused(checkpointname);
#endif
}
