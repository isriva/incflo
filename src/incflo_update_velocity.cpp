#include <incflo.H>

using namespace amrex;

void incflo::update_velocity (StepType step_type, Vector<MultiFab>& vel_eta, Vector<MultiFab>& vel_forces)
{
    BL_PROFILE("incflo::update_velocity");

    Real new_time = m_cur_time + m_dt;

    Real l_dt   = m_dt;
    Real l_half = Real(0.5);

    if (step_type == StepType::Predictor) {

        // *************************************************************************************
        // Define (or if advection_type != "MOL", re-define) the forcing terms, without the viscous terms
        //    and using the half-time density
        // *************************************************************************************
        compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_old_const(),
                           get_density_nph_const(), get_tracer_old_const(), get_tracer_new_const());
//        compute_stochastic_velocity_force(get_density_nph_const(), GetVecOfConstPtrs(vel_eta));
        add_stochastic_velocity_force(GetVecOfPtrs(vel_forces), step_type);

        // *************************************************************************************
        // Update the velocity
        // *************************************************************************************
        for (int lev = 0; lev <= finest_level; lev++)
        {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = ld.density_nph.const_array(mfi);

            if (m_diff_type == DiffusionType::Implicit) {

                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                    // Here divtau_o is the difference of tensor and scalar divtau_o!
                    if (m_advect_momentum) {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        });
                    }
                } else {
                    if (m_advect_momentum) {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                         vel(i,j,k,1) *= rho_old(i,j,k);,
                                         vel(i,j,k,2) *= rho_old(i,j,k););
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0));,
                                         vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1));,
                                         vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)););
                        });
                    }
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {

                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+l_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+l_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+l_half*divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+l_half*divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+l_half*divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+l_half*divtau_o(i,j,k,2)););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) *= rho_old(i,j,k);,
                                     vel(i,j,k,1) *= rho_old(i,j,k);,
                                     vel(i,j,k,2) *= rho_old(i,j,k););
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+rho_nph(i,j,k)*vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+rho_nph(i,j,k)*vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+rho_nph(i,j,k)*vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));,
                                     vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));,
                                     vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2)););
                    });
                }
            }
        } // mfi
        } // lev

    } else if (step_type == StepType::Corrector) {

        // *************************************************************************************
        // Define the forcing terms to use in the final update (using half-time density)
        // *************************************************************************************
        compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                           get_density_nph_const(), get_tracer_old_const(), get_tracer_new_const());
//        compute_stochastic_velocity_force(get_density_nph_const(), GetVecOfConstPtrs(vel_eta));
        add_stochastic_velocity_force(GetVecOfPtrs(vel_forces), step_type);

        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity.const_array(mfi);
            Array4<Real const> const& dvdt_o = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);

            Array4<Real const> const& rho_old  = ld.density_o.const_array(mfi);
            Array4<Real const> const& rho_new  = ld.density.const_array(mfi);
            Array4<Real const> const& rho_nph  = ld.density_nph.const_array(mfi);

            if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);

                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + l_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) =  rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                     l_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                   + l_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) =  rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                     l_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                   + l_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                   + rho_nph(i,j,k) * vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+  dvdt(i,j,k,0))
                                                  + l_half*(divtau_o(i,j,k,0)+divtau(i,j,k,0))
                                                  + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+  dvdt(i,j,k,1))
                                                  + l_half*(divtau_o(i,j,k,1)+divtau(i,j,k,1))
                                                  + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+  dvdt(i,j,k,2))
                                                  + l_half*(divtau_o(i,j,k,2)+divtau(i,j,k,2))
                                                  + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);

                if (m_advect_momentum) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                   +l_half*(divtau_o(i,j,k,0) ) + rho_nph(i,j,k)*vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                   +l_half*(divtau_o(i,j,k,1) ) + rho_nph(i,j,k)*vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                   +l_half*(divtau_o(i,j,k,2) ) + rho_nph(i,j,k)*vel_f(i,j,k,2) ););

                        AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                     vel(i,j,k,1) /= rho_new(i,j,k);,
                                     vel(i,j,k,2) /= rho_new(i,j,k););
                    });
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                  + l_half* divtau_o(i,j,k,0) + vel_f(i,j,k,0) );,
                                     vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                  + l_half* divtau_o(i,j,k,1) + vel_f(i,j,k,1) );,
                                     vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                    l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                  + l_half* divtau_o(i,j,k,2) + vel_f(i,j,k,2) ););
                    });
                }
            }
            else if (m_diff_type == DiffusionType::Implicit)
            {
                if (use_tensor_correction)
                {
                    Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);
                    if (m_advect_momentum)
                    {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) + divtau(i,j,k,0));,
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) + divtau(i,j,k,1));,
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) + divtau(i,j,k,2)););

                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0) +   dvdt(i,j,k,0))
                                                      + vel_f(i,j,k,0) + divtau(i,j,k,0));,
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(dvdt_o(i,j,k,1) +   dvdt(i,j,k,1))
                                                      + vel_f(i,j,k,1) + divtau(i,j,k,1) );,
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2) +   dvdt(i,j,k,2))
                                                      + vel_f(i,j,k,2) + divtau(i,j,k,2) ););
                        });
                    }
                } else {
                    if (m_advect_momentum)
                    {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = rho_old(i,j,k) * vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = rho_old(i,j,k) * vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = rho_old(i,j,k) * vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2))
                                                      + rho_nph(i,j,k) * vel_f(i,j,k,2) ););
                            AMREX_D_TERM(vel(i,j,k,0) /= rho_new(i,j,k);,
                                         vel(i,j,k,1) /= rho_new(i,j,k);,
                                         vel(i,j,k,2) /= rho_new(i,j,k););
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            AMREX_D_TERM(vel(i,j,k,0) = vel_o(i,j,k,0) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,0)+dvdt(i,j,k,0)) + vel_f(i,j,k,0) );,
                                         vel(i,j,k,1) = vel_o(i,j,k,1) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,1)+dvdt(i,j,k,1)) + vel_f(i,j,k,1) );,
                                         vel(i,j,k,2) = vel_o(i,j,k,2) + l_dt * (
                                                        l_half*(  dvdt_o(i,j,k,2)+dvdt(i,j,k,2)) + vel_f(i,j,k,2) ););
                        });
                    }
                }
            }
        }
        } // lev
    } else if (step_type == StepType::RK3StageOne || step_type == StepType::RK3StageTwo || step_type == StepType::RK3StageThree) {
        compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(), get_density_nph_const(),
                           get_tracer_old_const(), get_tracer_new_const(),
                           m_rk3_use_lagged_pressure_gradient);
        add_stochastic_velocity_force(GetVecOfPtrs(vel_forces), step_type);
        const bool stage_one = (step_type == StepType::RK3StageOne);
        const Real base = stage_one ? Real(1.0) :
                          (step_type == StepType::RK3StageTwo ? Real(0.75) : Real(1.0/3.0));
        const Real current = stage_one ? Real(1.0) :
                             (step_type == StepType::RK3StageTwo ? Real(0.25) : Real(2.0/3.0));
        const Real state_coeff = stage_one ? Real(0.0) : current;
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& vel = ld.velocity.array(mfi);
                Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
                Array4<Real const> const dvdt = stage_one ? ld.conv_velocity_o.const_array(mfi) : ld.conv_velocity.const_array(mfi);
                Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
                Array4<Real const> const divtau = stage_one ? ld.divtau_o.const_array(mfi) : ld.divtau.const_array(mfi);
                Array4<Real const> const& rho_nph = ld.density_nph.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    if (m_advect_momentum) {
                        AMREX_D_TERM(vel(i,j,k,0) = base*vel_o(i,j,k,0) + state_coeff*vel(i,j,k,0) + current*l_dt*((dvdt(i,j,k,0) + divtau(i,j,k,0))/rho_nph(i,j,k) + vel_f(i,j,k,0));,
                                     vel(i,j,k,1) = base*vel_o(i,j,k,1) + state_coeff*vel(i,j,k,1) + current*l_dt*((dvdt(i,j,k,1) + divtau(i,j,k,1))/rho_nph(i,j,k) + vel_f(i,j,k,1));,
                                     vel(i,j,k,2) = base*vel_o(i,j,k,2) + state_coeff*vel(i,j,k,2) + current*l_dt*((dvdt(i,j,k,2) + divtau(i,j,k,2))/rho_nph(i,j,k) + vel_f(i,j,k,2)););
                    } else {
                        AMREX_D_TERM(vel(i,j,k,0) = base*vel_o(i,j,k,0) + state_coeff*vel(i,j,k,0) + current*l_dt*(dvdt(i,j,k,0) + divtau(i,j,k,0) + vel_f(i,j,k,0));,
                                     vel(i,j,k,1) = base*vel_o(i,j,k,1) + state_coeff*vel(i,j,k,1) + current*l_dt*(dvdt(i,j,k,1) + divtau(i,j,k,1) + vel_f(i,j,k,1));,
                                     vel(i,j,k,2) = base*vel_o(i,j,k,2) + state_coeff*vel(i,j,k,2) + current*l_dt*(dvdt(i,j,k,2) + divtau(i,j,k,2) + vel_f(i,j,k,2)););
                    }
                });
            }
        }
    } // RK3 Stage

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            fillphysbc_density (lev, new_time, m_leveldata[lev]->density , ng_diffusion);
        }

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : l_half*m_dt;
        diffuse_velocity(get_velocity_new(), get_density_new(), GetVecOfConstPtrs(vel_eta), dt_diff);
    }
}

