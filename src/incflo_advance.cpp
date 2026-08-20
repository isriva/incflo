#include <incflo.H>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = static_cast<Real>(ParallelDescriptor::second());

    // Compute time step size
    int initialisation = ( m_dt < 0 );
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    // Set old/new times for fillpatching. Stage 1 produces U^(1), associated
    // with t^n + dt; Stage 2 produces U^(2), associated with t^n + dt/2.

    const bool rk3 = uses_RK3_timestepping();
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        m_t_old[lev] = m_cur_time;
        m_t_new[lev] = rk3
            ? m_cur_time + m_dt
            : m_cur_time + m_dt;
    }

    if (m_verbose > 0)
    {
        amrex::Print() << "\nStep " << m_nstep + 1
                       << ": from old_time " << m_cur_time
                       << " to new time " << m_cur_time + m_dt
                       << " with dt = " << m_dt << ".\n" << "\n";
    }

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();
    copy_from_new_to_old_temperature();

    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
        if (m_use_temperature) {
            fillpatch_temperature(lev, m_t_old[lev], m_leveldata[lev]->temperature_o, ng);
        }
    }

#ifdef AMREX_USE_EB
    if (m_eb_flow.enabled) {
       for (int lev = 0; lev <= finest_level; ++lev) {
         set_eb_velocity(lev, m_t_old[lev], *get_velocity_eb()[lev], 1);
         set_eb_density(lev, m_t_old[lev], *get_density_eb()[lev], 1);
       }
    }
    if (m_advect_tracer && !m_eb_flow.tracer.empty()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            set_eb_tracer(lev, m_t_old[lev], *get_tracer_eb()[lev], 1);
        }
    }
    if (m_use_temperature && !m_eb_flow.temperature.empty()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            set_eb_temperature(lev, m_t_old[lev], *get_temperature_eb()[lev], 1);
        }
    }
#endif

    ApplyPredictor(uses_RK3_timestepping() ? StepType::RK3StageOne : StepType::Predictor);

    if (uses_predictor_corrector_advection()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
            if (m_use_temperature) {
                fillpatch_temperature(lev, m_t_new[lev], m_leveldata[lev]->temperature, ng);
            }
        }

        ApplyCorrector(StepType::Corrector);
    }

    if (uses_RK3_timestepping()) {

        // Second stage
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
            if (m_use_temperature) {
                fillpatch_temperature(lev, m_t_new[lev], m_leveldata[lev]->temperature, ng);
            }
        }

        ApplyCorrector(StepType::RK3StageTwo);

        for (int lev = 0; lev <= finest_level; ++lev) {
            m_t_new[lev] = m_cur_time + m_dt / Real(2.0);
        }

        // Third stage
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
            if (m_use_temperature) {
                fillpatch_temperature(lev, m_t_new[lev], m_leveldata[lev]->temperature, ng);
            }
        }

        ApplyCorrector(StepType::RK3StageThree);
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_t_new[lev] = m_cur_time + m_dt;
        }
    }

    if (m_subgrid_modeling) {
        // Generate U'(t_n) using U'(t_{n-1}) and \bar{U}(t_n) 
        // Use this to approximate the turbulent stress tensor at time t_n
        compute_turb_stress_subgrid_force(m_subgrid_gamma_d, m_subgrid_gamma_f);
    }

#ifdef INCFLO_USE_PARTICLES
    particleData.Redistribute();
#endif

    // // This sums over all levels
    // if (m_test_tracer_conservation) {
    //     Real sum = volumeWeightedSum(get_tracer_new_const(),0,geom,ref_ratio);
    //     amrex::Print() << "Sum tracer volume wgt2 = " << m_cur_time+m_dt << " " <<
    //                        sum << "\n";
    // }

    // Stop timing current time step
    Real end_step = static_cast<Real>(ParallelDescriptor::second()) - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << "\n";
    }
}

