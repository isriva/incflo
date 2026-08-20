#include <incflo.H>

#include <AMReX_Random.H>

#include <cmath>

using namespace amrex;

void drift_diffusion_strain_ansatz(Vector<MultiFab const*> const& v_prime_displaced, Vector<MultiFab const*> const& v_bar_displaced, Real gamma_d)
{
    // Compute the large scale strain rate tensor S_{ij} = .5 (\partial_{x_j} \Bar{v}_i + \partial_{x_i} \Bar{v}_j)

    // Compute the strain rate magnitude \sqrt{2 S_{ij} S_{ij}}
}

void compute_turb_stress_subgrid_force (Real gamma_d, Real gamma_f, Vector<MultiFab const*> const& eta)
{
    if (!m_subgrid_modeling) {
        m_turb_stress_subgrid_force.clear();
        return;
    }

    BL_PROFILE("incflo::compute_turb_stress_subgrid_force");

    if (m_dt <= Real(0.0)) {
        amrex::Abort("turbulent stress tensor for subgrid modeling require a positive timestep");
    }

    m_turb_stress_subgrid_force.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        m_turb_stress_subgrid_force[lev].define(eta[lev]->boxArray(), eta[lev]->DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(),
                                           eta[lev]->Factory());
        m_turb_stress_subgrid_force[lev].setVal(0.0);

        auto& ld = *m_leveldata[lev];    
    
        // Interpolate velocity and subgrid_velocity to x* then apply MAC
        amrex::MultiFab subgrid_velocity_x_star(ld.subgrid_velocity.boxArray(), ld.subgrid_velocity.DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(), ld.subgrid_velocity.Factory());
        amrex::MultiFab velocity_x_star(ld.velocity.boxArray(), ld.velocity.DistributionMap(),
                                           AMREX_SPACEDIM, nghost_force(), MFInfo(), ld.velocity.Factory());

        // Compute drift term a(v', \bar{v}, gamma_d) and diffusion term b(v', \bar{v}, gamma_f)
        if (m_subgrid_model_type == "strain_ansatz"){
            drift_diffusion_strain_ansatz(Vector<MultiFab const*> const& v_prime_displaced, Vector<MultiFab const*> const& v_bar_displaced, Real gamma_d)
        }

        // Sample dW

        // Compute dv' = a(v', \bar{v}, gamma_d) dt + b(v', \bar{v}, gamma_f) dW

        // Form v'(x,t) = v'(x^*, t - dt) + dx'

        // Form \tilde{v} = \bar{v} + v'

        // Compute tau = \bar{\tilde{v} \tilde{v}} - \bar{v} \bar{v}

    }

}

void add_turb_stress_subgrid_force (Vector<MultiFab*> const& vel_forces, StepType step_type) const
{
    if (!m_subgrid_modeling) {
        return;
    }

    AMREX_ALWAYS_ASSERT(m_turb_stress_subgrid_force.size() == vel_forces.size());

    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Add(*vel_forces[lev], m_turb_stress_subgrid_force[lev],
                      0, 0, AMREX_SPACEDIM, 0);
    }
}
