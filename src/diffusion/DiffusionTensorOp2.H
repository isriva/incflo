#ifndef DIFFUSION_TENSOR_OP2_H_
#define DIFFUSION_TENSOR_OP2_H_

#include <AMReX_MLMG.H>

#include <MyTensorOp.H>

//
// Solver for the implicit part of the diffusion equation:
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//

class incflo;

class DiffusionTensorOp2
{
public:
    DiffusionTensorOp2 (incflo* a_incflo);

    void diffuse_velocity (amrex::Vector<amrex::MultiFab*> const& velocity,
                           amrex::Vector<amrex::MultiFab*> const& density,
                           amrex::Vector<amrex::MultiFab const*> const& eta,
                           amrex::Real dt);

    void compute_divtau (amrex::Vector<amrex::MultiFab*> const& divtau,
                         amrex::Vector<amrex::MultiFab const*> const& velocity,
                         amrex::Vector<amrex::MultiFab const*> const& density,
                         amrex::Vector<amrex::MultiFab const*> const& eta);

private:

    void readParameters ();

    incflo* m_incflo;

    std::unique_ptr<amrex::MyTensorOp> m_reg_solve_op2;
    std::unique_ptr<amrex::MyTensorOp> m_reg_apply_op2;

    // DiffusionOp verbosity
    int m_verbose = 0;

    // Options to control MLMG behavior
    int m_mg_verbose = 0;
    int m_mg_bottom_verbose = 0;
    int m_mg_max_iter = 100;
    int m_mg_bottom_maxiter = 100;
    int m_mg_max_fmg_iter = 0;
    int m_mg_max_coarsening_level = 100;
    int m_mg_maxorder = 2;
    int m_num_pre_smooth = 2;
    int m_num_post_smooth = 2;
#ifdef AMREX_USE_FLOAT
    amrex::Real m_mg_rtol = 1.0e-4;
    amrex::Real m_mg_atol = 1.0e-7;
#else
    amrex::Real m_mg_rtol = 1.0e-11;
    amrex::Real m_mg_atol = 1.0e-14;
#endif
    std::string m_bottom_solver = "bicgstab";
};

#endif
