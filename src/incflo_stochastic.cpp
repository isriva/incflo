#include <incflo.H>

#include <AMReX_Random.H>

#include <cmath>

using namespace amrex;

void
incflo::compute_stochastic_velocity_force (Vector<MultiFab const*> const& density,
                                           Vector<MultiFab const*> const& eta)
{
    if (!m_use_stochastic_velocity_fluxes) {
        m_stochastic_vel_force.clear();
        return;
    }

    BL_PROFILE("incflo::compute_stochastic_velocity_force");
    amrex::ignore_unused(density);

    if (m_dt <= Real(0.0)) {
        amrex::Abort("stochastic velocity fluxes require a positive timestep");
    }

    m_stochastic_vel_force.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        m_stochastic_vel_force[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(),
                                           eta[lev]->Factory());
        m_stochastic_vel_force[lev].setVal(0.0);

        Real cell_volume = Real(1.0);
        auto const dx = geom[lev].CellSizeArray();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            cell_volume *= dx[dir];
        }
#if (AMREX_SPACEDIM == 2)
        cell_volume *= m_stochastic_cell_depth;
#endif

        Array<MultiFab, AMREX_SPACEDIM> eta_face = average_velocity_eta_to_faces(lev, *eta[lev]);
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux;

        const Real coeff = std::sqrt(Real(2.0) * m_stochastic_k_B *
                                     m_stochastic_temperature / (cell_volume * m_dt));

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            stochastic_flux[dir].define(eta_face[dir].boxArray(), eta_face[dir].DistributionMap(),
                                        AMREX_SPACEDIM, 0, MFInfo(), eta_face[dir].Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(stochastic_flux[dir], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& flux = stochastic_flux[dir].array(mfi);
                Array4<Real const> const& face_eta = eta_face[dir].const_array(mfi);

                ParallelForRNG(bx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                       RandomEngine const& engine) noexcept
                {
                    flux(i,j,k,n) = coeff * std::sqrt(face_eta(i,j,k)) *
                                    RandomNormal(Real(0.0), Real(1.0), engine);
                });
            }

            stochastic_flux[dir].OverrideSync(geom[lev].periodicity());
            
            // stochastic_flux[dir].EnforcePeriodicity(0, AMREX_SPACEDIM,
            //                                         geom[lev].periodicity());
        }

        Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux_ptrs[dir] = &stochastic_flux[dir];
        }
        amrex::computeDivergence(m_stochastic_vel_force[lev], flux_ptrs, geom[lev]);

        const Real rhoinv = Real(1.0) / m_ro_0;
        m_stochastic_vel_force[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
        m_stochastic_vel_force[lev].FillBoundary(geom[lev].periodicity());
    }
}

void
incflo::add_stochastic_velocity_force (Vector<MultiFab*> const& vel_forces) const
{
    if (!m_use_stochastic_velocity_fluxes) {
        return;
    }

    AMREX_ALWAYS_ASSERT(m_stochastic_vel_force.size() == vel_forces.size());

    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Add(*vel_forces[lev], m_stochastic_vel_force[lev],
                      0, 0, AMREX_SPACEDIM, 0);
    }
}
