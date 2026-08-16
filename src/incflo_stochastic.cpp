#include <incflo.H>

#include <AMReX_Random.H>

#include <cmath>

using namespace amrex;

void
incflo::compute_stochastic_velocity_force (Vector<MultiFab const*> const& density,
                                           Vector<MultiFab const*> const& eta)
{
    if (!m_use_stochastic_velocity_fluxes) {
        m_stochastic_vel_force_pred.clear();
        m_stochastic_vel_force_corr.clear();
        return;
    }

    BL_PROFILE("incflo::compute_stochastic_velocity_force");
    amrex::ignore_unused(density);

    if (m_dt <= Real(0.0)) {
        amrex::Abort("stochastic velocity fluxes require a positive timestep");
    }

    m_stochastic_vel_force_pred.resize(finest_level+1);
    m_stochastic_vel_force_corr.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        m_stochastic_vel_force_pred[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(),
                                           eta[lev]->Factory());
        m_stochastic_vel_force_pred[lev].setVal(0.0);
        m_stochastic_vel_force_corr[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(),
                                           eta[lev]->Factory());
        m_stochastic_vel_force_corr[lev].setVal(0.0);

        Real cell_volume = Real(1.0);
        auto const dx = geom[lev].CellSizeArray();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            cell_volume *= dx[dir];
        }
#if (AMREX_SPACEDIM == 2)
        cell_volume *= m_stochastic_cell_depth;
#endif

        Array<MultiFab, AMREX_SPACEDIM> eta_face = average_velocity_eta_to_faces(lev, *eta[lev]);
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_A;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_B;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_pred;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_corr;

        const Real coeff = std::sqrt(Real(2.0) * m_stochastic_k_B *
                                     m_stochastic_temperature / (cell_volume * m_dt));

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            stochastic_flux_A[dir].define(eta_face[dir].boxArray(), eta_face[dir].DistributionMap(),
                                        AMREX_SPACEDIM, 0, MFInfo(), eta_face[dir].Factory());
            stochastic_flux_B[dir].define(eta_face[dir].boxArray(), eta_face[dir].DistributionMap(),
                                        AMREX_SPACEDIM, 0, MFInfo(), eta_face[dir].Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(stochastic_flux_A[dir], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& flux_A = stochastic_flux_A[dir].array(mfi);
                Array4<Real> const& flux_B = stochastic_flux_B[dir].array(mfi);
                Array4<Real const> const& face_eta = eta_face[dir].const_array(mfi);

                ParallelForRNG(bx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                       RandomEngine const& engine) noexcept
                {
                    flux_A(i,j,k,n) = coeff * std::sqrt(face_eta(i,j,k)) *
                                    RandomNormal(Real(0.0), Real(1.0), engine);
                    flux_B(i,j,k,n) = coeff * std::sqrt(face_eta(i,j,k)) *
                                    RandomNormal(Real(0.0), Real(1.0), engine);
                });
            }

            stochastic_flux_A[dir].OverrideSync(geom[lev].periodicity());
            stochastic_flux_B[dir].OverrideSync(geom[lev].periodicity());
       }

       // Predictor
       {
           Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
           for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
               stochastic_flux_pred[dir].define(stochastic_flux_A[dir].boxArray(),
                       stochastic_flux_A[dir].DistributionMap(),
                       AMREX_SPACEDIM, 0, MFInfo(), stochastic_flux_A[dir].Factory());
               MultiFab::LinComb(stochastic_flux_pred[dir],
                       m_stoch_weights[0], stochastic_flux_A[dir], 0,
                       m_stoch_weights[1], stochastic_flux_B[dir], 0,
                       0, AMREX_SPACEDIM, 0);
               flux_ptrs[dir] = &stochastic_flux_pred[dir];
           }
           amrex::computeDivergence(m_stochastic_vel_force_pred[lev], flux_ptrs, geom[lev]);

           const Real rhoinv = Real(1.0) / m_ro_0;
           m_stochastic_vel_force_pred[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
           m_stochastic_vel_force_pred[lev].FillBoundary(geom[lev].periodicity());
       }

       // Corrector
       {
           Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
           for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
               stochastic_flux_corr[dir].define(stochastic_flux_A[dir].boxArray(),
                       stochastic_flux_A[dir].DistributionMap(),
                       AMREX_SPACEDIM, 0, MFInfo(), stochastic_flux_A[dir].Factory());
               MultiFab::LinComb(stochastic_flux_corr[dir],
                       m_stoch_weights[2], stochastic_flux_A[dir], 0,
                       m_stoch_weights[3], stochastic_flux_B[dir], 0,
                       0, AMREX_SPACEDIM, 0);
               flux_ptrs[dir] = &stochastic_flux_corr[dir];
           }
           amrex::computeDivergence(m_stochastic_vel_force_corr[lev], flux_ptrs, geom[lev]);

           const Real rhoinv = Real(1.0) / m_ro_0;
           m_stochastic_vel_force_corr[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
           m_stochastic_vel_force_corr[lev].FillBoundary(geom[lev].periodicity());
       }
    }
}

void
incflo::compute_stochastic_velocity_force_RK3 (Vector<MultiFab const*> const& density,
                                               Vector<MultiFab const*> const& eta)
{
    if (!m_use_stochastic_velocity_fluxes) {
        m_stochastic_vel_force_RK3_stage_one.clear();
        m_stochastic_vel_force_RK3_stage_two.clear();
        m_stochastic_vel_force_RK3_stage_three.clear();
        return;
    }

    BL_PROFILE("incflo::compute_stochastic_velocity_force_RK3");
    amrex::ignore_unused(density);
    if (m_dt <= Real(0.0)) {
        amrex::Abort("stochastic velocity fluxes require a positive timestep");
    }

    const Real beta[3] = {
        (Real(2.0)*std::sqrt(Real(2.0)) + std::sqrt(Real(3.0))) / Real(5.0),
        (-Real(4.0)*std::sqrt(Real(2.0)) + Real(3.0)*std::sqrt(Real(3.0))) / Real(5.0),
        (std::sqrt(Real(2.0)) - Real(2.0)*std::sqrt(Real(3.0))) / Real(10.0)
    };
    m_stochastic_vel_force_RK3_stage_one.resize(finest_level+1);
    m_stochastic_vel_force_RK3_stage_two.resize(finest_level+1);
    m_stochastic_vel_force_RK3_stage_three.resize(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        Real cell_volume = Real(1.0);
        auto const dx = geom[lev].CellSizeArray();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) cell_volume *= dx[dir];
#if (AMREX_SPACEDIM == 2)
        cell_volume *= m_stochastic_cell_depth;
#endif

        m_stochastic_vel_force_RK3_stage_one[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                                       AMREX_SPACEDIM, nghost_force(), MFInfo(), eta[lev]->Factory());
        m_stochastic_vel_force_RK3_stage_one[lev].setVal(0.0);
        m_stochastic_vel_force_RK3_stage_two[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                                       AMREX_SPACEDIM, nghost_force(), MFInfo(), eta[lev]->Factory());
        m_stochastic_vel_force_RK3_stage_two[lev].setVal(0.0);
        m_stochastic_vel_force_RK3_stage_three[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                                       AMREX_SPACEDIM, nghost_force(), MFInfo(), eta[lev]->Factory());
        m_stochastic_vel_force_RK3_stage_three[lev].setVal(0.0);

        Array<MultiFab, AMREX_SPACEDIM> eta_face = average_velocity_eta_to_faces(lev, *eta[lev]);
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_A;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_B;
        const Real coeff = std::sqrt(Real(2.0) * m_stochastic_k_B *
                                     m_stochastic_temperature / (cell_volume * m_dt));

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            stochastic_flux_A[dir].define(eta_face[dir].boxArray(), eta_face[dir].DistributionMap(),
                               AMREX_SPACEDIM, 0, MFInfo(), eta_face[dir].Factory());
            stochastic_flux_B[dir].define(eta_face[dir].boxArray(), eta_face[dir].DistributionMap(),
                               AMREX_SPACEDIM, 0, MFInfo(), eta_face[dir].Factory());
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(stochastic_flux_A[dir], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& flux_A = stochastic_flux_A[dir].array(mfi);
                Array4<Real> const& flux_B = stochastic_flux_B[dir].array(mfi);
                Array4<Real const> const& face_eta = eta_face[dir].const_array(mfi);
                ParallelForRNG(bx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                      RandomEngine const& engine) noexcept {
                    flux_A(i,j,k,n) = coeff * std::sqrt(face_eta(i,j,k)) * RandomNormal(Real(0.0), Real(1.0), engine);
                    flux_B(i,j,k,n) = coeff * std::sqrt(face_eta(i,j,k)) * RandomNormal(Real(0.0), Real(1.0), engine);
                });
            }
            stochastic_flux_A[dir].OverrideSync(geom[lev].periodicity());
            stochastic_flux_B[dir].OverrideSync(geom[lev].periodicity());
        }

        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_RK3_one;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_RK3_two;
        Array<MultiFab, AMREX_SPACEDIM> stochastic_flux_RK3_three;

        // RK3 stage one
        {
            Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                stochastic_flux_RK3_one[dir].define(stochastic_flux_A[dir].boxArray(),
                        stochastic_flux_A[dir].DistributionMap(),
                        AMREX_SPACEDIM, 0, MFInfo(), stochastic_flux_A[dir].Factory());
                MultiFab::LinComb(stochastic_flux_RK3_one[dir], Real(1.0), stochastic_flux_A[dir], 0,
                                  beta[0], stochastic_flux_B[dir], 0, 0, AMREX_SPACEDIM, 0);
                flux_ptrs[dir] = &stochastic_flux_RK3_one[dir];
            }
            amrex::computeDivergence(m_stochastic_vel_force_RK3_stage_one[lev], flux_ptrs, geom[lev]);

            const Real rhoinv = Real(1.0) / m_ro_0;
            m_stochastic_vel_force_RK3_stage_one[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
            m_stochastic_vel_force_RK3_stage_one[lev].FillBoundary(geom[lev].periodicity());
        }

        // RK3 stage two
        {
            Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                stochastic_flux_RK3_two[dir].define(stochastic_flux_A[dir].boxArray(),
                        stochastic_flux_A[dir].DistributionMap(),
                        AMREX_SPACEDIM, 0, MFInfo(), stochastic_flux_A[dir].Factory());
                MultiFab::LinComb(stochastic_flux_RK3_two[dir], Real(1.0), stochastic_flux_A[dir], 0,
                                  beta[1], stochastic_flux_B[dir], 0, 0, AMREX_SPACEDIM, 0);
                flux_ptrs[dir] = &stochastic_flux_RK3_two[dir];
            }
            amrex::computeDivergence(m_stochastic_vel_force_RK3_stage_two[lev], flux_ptrs, geom[lev]);
            const Real rhoinv = Real(1.0) / m_ro_0;
            m_stochastic_vel_force_RK3_stage_two[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
            m_stochastic_vel_force_RK3_stage_two[lev].FillBoundary(geom[lev].periodicity());
        }

        // RK3 stage three
        {
            Array<MultiFab const*, AMREX_SPACEDIM> flux_ptrs;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                stochastic_flux_RK3_three[dir].define(stochastic_flux_A[dir].boxArray(),
                        stochastic_flux_A[dir].DistributionMap(),
                        AMREX_SPACEDIM, 0, MFInfo(), stochastic_flux_A[dir].Factory());
                MultiFab::LinComb(stochastic_flux_RK3_three[dir], Real(1.0), stochastic_flux_A[dir], 0,
                                  beta[2], stochastic_flux_B[dir], 0, 0, AMREX_SPACEDIM, 0);
                flux_ptrs[dir] = &stochastic_flux_RK3_three[dir];
            }
            amrex::computeDivergence(m_stochastic_vel_force_RK3_stage_three[lev], flux_ptrs, geom[lev]);
            const Real rhoinv = Real(1.0) / m_ro_0;
            m_stochastic_vel_force_RK3_stage_three[lev].mult(rhoinv, 0, AMREX_SPACEDIM, 0);
            m_stochastic_vel_force_RK3_stage_three[lev].FillBoundary(geom[lev].periodicity());
        }
    }
}

void
incflo::add_stochastic_velocity_force (Vector<MultiFab*> const& vel_forces, StepType step_type) const
{
    if (!m_use_stochastic_velocity_fluxes) {
        return;
    }

    if ((step_type == StepType::Predictor && !uses_RK3_timestepping()) ||
        step_type == StepType::Corrector) {
        AMREX_ALWAYS_ASSERT(m_stochastic_vel_force_pred.size() == vel_forces.size());
        AMREX_ALWAYS_ASSERT(m_stochastic_vel_force_corr.size() == vel_forces.size());
    } else {
        AMREX_ALWAYS_ASSERT(m_stochastic_vel_force_RK3_stage_one.size() == vel_forces.size());
        AMREX_ALWAYS_ASSERT(m_stochastic_vel_force_RK3_stage_two.size() == vel_forces.size());
        AMREX_ALWAYS_ASSERT(m_stochastic_vel_force_RK3_stage_three.size() == vel_forces.size());
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        if (step_type == StepType::Predictor) {
            if (uses_RK3_timestepping()) {
                MultiFab::Add(*vel_forces[lev], m_stochastic_vel_force_RK3_stage_one[lev],
                              0, 0, AMREX_SPACEDIM, 0);
                continue;
            }
            MultiFab::Add(*vel_forces[lev], m_stochastic_vel_force_pred[lev],
                          0, 0, AMREX_SPACEDIM, 0);
        } else if (step_type == StepType::Corrector) {
            MultiFab::Add(*vel_forces[lev], m_stochastic_vel_force_corr[lev],
                          0, 0, AMREX_SPACEDIM, 0);
        } else if (step_type == StepType::RK3StageOne ||
                   step_type == StepType::RK3StageTwo ||
                   step_type == StepType::RK3StageThree) {
            int stage = static_cast<int>(step_type) -
                        static_cast<int>(StepType::RK3StageOne);
            Vector<MultiFab> const* stage_force =
                (stage == 0) ? &m_stochastic_vel_force_RK3_stage_one :
                (stage == 1) ? &m_stochastic_vel_force_RK3_stage_two :
                               &m_stochastic_vel_force_RK3_stage_three;
            MultiFab::Add(*vel_forces[lev], (*stage_force)[lev],
                          0, 0, AMREX_SPACEDIM, 0);
        }
    }
}
