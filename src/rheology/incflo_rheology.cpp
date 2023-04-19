#include <incflo.H>
#include <incflo_derive_K.H>

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real expterm (amrex::Real nu) noexcept
{
    return (nu < Real(1.e-9)) ? (Real(1.0)-Real(0.5)*nu+nu*nu*Real(1.0/6.0)-(nu*nu*nu)*Real(1./24.))
                        : -std::expm1(-nu)/nu;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertialNum (amrex::Real sr, amrex::Real p_ext, amrex::Real ro_0, amrex::Real diam, amrex::Real mu_1, amrex::Real A_1, amrex::Real alpha_1) noexcept
{
    return mu_1 + A_1*std::pow((sr/2)*diam/(std::pow(p_ext/ro_0, 0.5)), alpha_1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_concentration (amrex::Real dens, amrex::Real rho0, amrex::Real rho1) noexcept
{
    return ((rho1/dens) - 1.0)/((rho1/rho0) - 1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixture_viscosity (amrex::Real conc, amrex::Real visc0, amrex::Real visc1) noexcept
{
    return 1.0/((conc/visc0) + ((1.0-conc)/visc1));
    //return visc1/(1.0 + conc*((visc1/visc0) - 1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::tuple<amrex::Real, bool> Viscosity_Single(const amrex::Real sr, const int order, const incflo::FLUID_t& fluid, const amrex::Real hyd_press, const amrex::Real gravity) 
{
    amrex::Real visc = 0.0;
    bool include = false;
    
    if (fluid.fluid_model == incflo::FluidModel::Newtonian) {
        if (order == 0) {visc = fluid.mu; include = true;}
        else {visc = 0.0; include = false;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::Powerlaw) {
        if (order == 0) {visc = (fluid.mu * std::pow(sr,fluid.n_0-1.0)); include = true;}
        else if (order == 1) {visc = (fluid.mu_1 * std::pow(sr,fluid.n_1-1.0)); include = true;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::Bingham) {
        if (order == 0) {visc = fluid.mu + fluid.tau_0 * expterm(sr/fluid.papa_reg) / fluid.papa_reg; include = true;}
        else if (order == 1) {visc = fluid.mu_1 + fluid.tau_1 * expterm(sr/fluid.papa_reg_1) / fluid.papa_reg_1; include = true;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::HerschelBulkley) {
        if (order == 0) {
            visc = (fluid.mu*std::pow(sr,fluid.n_0)+fluid.tau_0)*expterm(sr/fluid.papa_reg)/fluid.papa_reg;
            include = true;
        }
        else if (order == 1) {
            visc = (fluid.mu_1*std::pow(sr,fluid.n_1)+fluid.tau_1)*expterm(sr/fluid.papa_reg_1)/fluid.papa_reg_1;
            include = true;
        }
    }
    else if (fluid.fluid_model == incflo::FluidModel::deSouzaMendesDutra) {
        if (order == 0) {
            visc = (fluid.mu*std::pow(sr,fluid.n_0)+fluid.tau_0)*expterm(sr*(fluid.eta_0/fluid.tau_0))*(fluid.eta_0/fluid.tau_0);
            include = true;
        }
        else if (order == 1) {
            visc = (fluid.mu_1*std::pow(sr,fluid.n_1)+fluid.tau_1)*expterm(sr*(fluid.eta_1/fluid.tau_1))*(fluid.eta_1/fluid.tau_1);
            include = true;
        }
    }
    else if (fluid.fluid_model == incflo::FluidModel::Granular) {
        if (order == 0) {
            amrex::Real hyd_press_reg = 0.5*(hyd_press + sqrt(hyd_press*hyd_press + fluid.papa_reg_press*fluid.papa_reg_press)); // regularized pressure (always positive)
            amrex::Real min_visc = Real(0.0);
            amrex::Real compute_visc = (fluid.tau_0*hyd_press_reg + fluid.A_0*std::pow(fluid.diam,fluid.alpha_0)*std::pow(fluid.rho,0.5*fluid.alpha_0)*std::pow(hyd_press_reg,1.0-0.5*fluid.alpha_0)*std::pow(sr,fluid.alpha_0))*expterm(sr/fluid.papa_reg)/fluid.papa_reg;
            visc = std::max(min_visc, compute_visc);
            include = true;
        }
        else if (order == 1) {
            visc = std::pow(2*(expterm(sr/fluid.papa_reg_1) / fluid.papa_reg_1),2)*(fluid.papa_reg_press)*inertialNum(sr,fluid.papa_reg_press, fluid.rho, fluid.diam, fluid.mu_1, fluid.A_1, 2*fluid.alpha_1);

            {
            amrex::Print() << "viscosity order 1:  " << visc << ", "  << std::endl;
            }
        
        }
        //else if (order == 2) {
        //    visc = -1*std::pow(2*(expterm(sr/papa_reg) / papa_reg),2)*(p_bg)*inertialNum(sr, p_bg, ro_0, diam, mu_3, A_3, 2*alpha_3);
        //}
    }
    return {visc, include};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Viscosity_VOF(const amrex::Real sr, const amrex::Real dens, const int order, const amrex::Vector<incflo::FLUID_t>& fluid_vof, const amrex::Real hyd_press, const amrex::Real gravity) 
{
    amrex::Real conc = get_concentration(dens,fluid_vof[0].rho,fluid_vof[1].rho);
    auto [visc0, include0] = Viscosity_Single(sr,order,fluid_vof[0],hyd_press,gravity);
    auto [visc1, include1] = Viscosity_Single(sr,order,fluid_vof[1],hyd_press,gravity);

    amrex::Real visc;
    if (include0 and include1) { // both fluid0 and fluid1 have rheologies of this order
        visc = mixture_viscosity(conc,visc0,visc1);
    }
    else if (include0 and !include1) { // only fluid0 rheology is  used for this order
        if (conc < 1e-12) visc = 0.0;
        else if (conc > 1.0-1e-12) visc = visc0;
        else visc = visc0*conc;
    }
    else if (!include0 and include1) { // only fluid1 rheology is used for this order
        if (conc < 1e-12) visc = visc1;
        else if (conc > 1.0-1e-12) visc = 0.0;
        else visc = visc1*(1.0 - conc);
    }
    else if (!include0 and !include1) { // neither fluid0 & fluid1 rheology are used for this order
        visc = 0.0;
    }

    return visc;

}

}

void incflo::compute_viscosity (amrex::Vector<amrex::MultiFab*> const& vel_eta,
                                amrex::Vector<amrex::MultiFab*> const& rho,
                                amrex::Vector<amrex::MultiFab*> const& vel,
                                amrex::Vector<amrex::MultiFab const*> const& p_nodal,
                                amrex::Real time, int nghost, int order)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        compute_viscosity_at_level(lev, vel_eta[lev], rho[lev], vel[lev], p_nodal[lev], geom[lev], time, nghost, order);
    }
}

#ifdef AMREX_USE_EB
void incflo::compute_viscosity_at_level (int lev,
#else
void incflo::compute_viscosity_at_level (int /*lev*/,
#endif
                                         amrex::MultiFab* vel_eta,
                                         amrex::MultiFab* rho,
                                         amrex::MultiFab* vel,
                                         amrex::MultiFab const* p_nodal,
                                         amrex::Geometry& lev_geom,
                                         amrex::Real /*time*/, int /*nghost*/, int order)
{

#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

    auto const& dx     = lev_geom.CellSizeArray();
    auto const& problo = lev_geom.ProbLoArray();
    auto const& probhi = lev_geom.ProbHiArray();
    Real idx           = Real(1.0) / lev_geom.CellSize(0);
    Real idy           = Real(1.0) / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
    Real idz           = Real(1.0) / lev_geom.CellSize(2);
#endif
        const Box& domain = lev_geom.Domain();
        const Dim3 domlo = amrex::lbound(domain);
        const Dim3 domhi = amrex::ubound(domain);
        Vector<Array<int,2>> bc_type(AMREX_SPACEDIM);
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = m_bc_type[ori];
            if (bct == BC::no_slip_wall) {
                if (side == Orientation::low)  bc_type[dir][0] = 2; 
                if (side == Orientation::high) bc_type[dir][1] = 2; 
            }
            else if (bct == BC::slip_wall) {
                if (side == Orientation::low)  bc_type[dir][0] = 1; 
                if (side == Orientation::high) bc_type[dir][1] = 1; 
            }
            else {
                if (side == Orientation::low)  bc_type[dir][0] = 0; 
                if (side == Orientation::high) bc_type[dir][1] = 0; 
            }
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& eta_arr = vel_eta->array(mfi);
        Array4<Real const> const& vel_arr = vel->const_array(mfi);
        Array4<Real const> const& rho_arr = rho->const_array(mfi);
        Array4<Real const> const& p_nodal_arr = p_nodal->const_array(mfi);
#ifdef AMREX_USE_EB
        auto const& flag_fab = flags[mfi];
        auto typ = flag_fab.getType(bx);
        if (typ == FabType::covered)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                eta_arr(i,j,k) = Real(0.0);
            });
        }
        else if (typ == FabType::singlevalued)
        {
            auto const& flag_arr = flag_fab.const_array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real sr = incflo_strainrate_eb(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr,flag_arr(i,j,k));
                Real dens = rho_arr(i,j,k);
                if (m_do_vof) {
                    eta_arr(i,j,k) = Viscosity_VOF(sr,dens,order,m_fluid_vof);
                }
                else {
                    auto [visc, include] = Viscosity_Single(sr,order,m_fluid);    
                    if (include) eta_arr(i,j,k) = visc;
                    else eta_arr(i,j,k) = 0.0;
                }
            });
        }
        else
#endif
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real sr = incflo_strainrate_nodal(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr,domlo,domhi,bc_type); //  get nodal strainrate
                Real dens = incflo_rho_nodal(i,j,k,rho_arr); // get nodal density
                Real depth_from_surface = probhi[2] - Real(k)*dx[2]; // get depth from the surface 
                Real hyd_press = m_p_amb_surface - m_gravity[2]*dens*depth_from_surface; // get hydrostatic pressure
                if (m_include_perturb_pressure) {
                   hyd_press += p_nodal_arr(i,j,k);                    
                }
                // nodal viscosity
                if (m_do_vof) {
                    eta_arr(i,j,k) = Viscosity_VOF(sr,dens,order,m_fluid_vof,hyd_press,m_gravity[2]);
                }
                else {
                    auto [visc, include] = Viscosity_Single(sr,order,m_fluid,hyd_press,m_gravity[2]);
                    if (include) eta_arr(i,j,k) = visc;
                    else eta_arr(i,j,k) = 0.0;
                }
            });
        }
    }
}

void incflo::compute_tracer_diff_coeff (Vector<MultiFab*> const& tra_eta, int nghost)
{
    for (auto *mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);
        }
    }
}
