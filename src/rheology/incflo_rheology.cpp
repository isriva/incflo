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
amrex::Real get_concentration (amrex::Real dens, amrex::Real rho0, amrex::Real rho1) noexcept
{
    return ((rho1/dens) - 1.0)/((rho1/rho0) - 1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixture_viscosity (amrex::Real conc, amrex::Real visc0, amrex::Real visc1) noexcept
{
    return visc1/(1.0 + conc*((visc1/visc0) - 1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertialNum (amrex::Real sr, amrex::Real prs, amrex::Real ro, amrex::Real diam, amrex::Real mu, amrex::Real A, amrex::Real alpha) noexcept
{
    return mu + A*std::pow((sr/2)*diam/(std::pow(prs/ro, 0.5)), alpha);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real eta1 (amrex::Real A, amrex::Real d, amrex::Real ro, amrex::Real alpha, amrex::Real sr, amrex::Real papa_reg, amrex::Real p)
{
    return (A*std::pow(d * std::pow(ro,1/2), alpha)) * std::pow(2*(expterm(sr/papa_reg)/papa_reg),-alpha+1) * std::pow(p, 1-alpha/2);
} 

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real eta2 (amrex::Real A, amrex::Real d, amrex::Real ro, amrex::Real alpha, amrex::Real sr, amrex::Real papa_reg, amrex::Real p)
{
    return (A*std::pow(d * std::pow(ro,1/2), alpha)) * std::pow(2*(expterm(sr/papa_reg)/papa_reg),-alpha+2) * std::pow(p, 1-alpha/2);
} 

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kappaterm (amrex::Real mu, amrex::Real p)
{
    return mu * p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::tuple<amrex::Real, bool> Viscosity_Single(const amrex::Real sr, const int order, const incflo::FLUID_t& fluid) 
{
    amrex::Real visc = 0.0;
    bool include = false;
    
    if (fluid.fluid_model == incflo::FluidModel::Newtonian) {
        if (order == 0) {visc = fluid.mu; include = true;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::Powerlaw) {
        if (order == 0) {
            if (sr < 1.e-15) {visc = fluid.mu; include = true;}
            else {visc = (fluid.mu * std::pow(sr,fluid.n_0-1.0)); include = true;}
        }
        else if (order == 1) {
            if (sr < 1.e-15) {visc = fluid.mu_1; include = true;}
            else {visc = (fluid.mu_1 * std::pow(sr,fluid.n_1-1.0)); include = true;}
        }
    }
    else if (fluid.fluid_model == incflo::FluidModel::Bingham) {
        if (order == 0) {visc = fluid.mu + fluid.tau_0 * expterm(sr/fluid.papa_reg) / fluid.papa_reg; include = true;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::HerschelBulkley) {
        if (order == 0) {
            //visc = (fluid.mu*std::pow(sr,fluid.n_0)+fluid.tau_0)*expterm(sr/fluid.papa_reg)/fluid.papa_reg;
            // return (mu*std::pow(sr,n_flow)+tau_0)*expterm(sr/papa_reg)/papa_reg;
            visc = ( fluid.mu*std::pow(sr,fluid.n_0-1.0) + (fluid.tau_0/sr)*(1.0-expterm(-1.0*sr/fluid.papa_reg)));
        }
        else if (order == 1) {
            // visc =  (mu_1*std::pow(sr,n_flow_1)+tau_1)*expterm(sr/papa_reg_1)/papa_reg_1;
            visc = ( fluid.mu_1*std::pow(sr,fluid.n_1-1.0) + (fluid.tau_1/sr)*(1.0-expterm(-1.0*sr/fluid.papa_reg_1)));
        }
    }
    else if (fluid.fluid_model == incflo::FluidModel::deSouzaMendesDutra) {
        visc = (fluid.mu*std::pow(sr,fluid.n_0)+fluid.tau_0)*expterm(sr*(fluid.eta_0/fluid.tau_0))*(fluid.eta_0/fluid.tau_0);
    }
    else if (fluid.fluid_model == incflo::FluidModel::Granular) {
        // If you want a pressure gradient due to gravity, add p_ext to p_bg
        // For the strainrate, power is zero for an initial test. Make it "1" for a real simulation
        Real p_bg = 1.0;
        if (order == 0) {
            visc = std::pow(2*(expterm(sr/fluid.papa_reg) / fluid.papa_reg),1)*(p_bg)*inertialNum(sr,p_bg, fluid.rho, fluid.diam, fluid.mu, fluid.A_0, fluid.alpha_0);
        }
        //else if (order == 1) {
        //    visc = std::pow(2*(expterm(sr/fluid.papa_reg_1) / fluid.papa_reg_1),2)*(fluid.p_bg)*inertialNum(sr,fluid.p_bg, fluid.rho, fluid.diam, fluid.mu_2, fluid.A_2, 2*fluid.alpha_2);
        //}
        //else if (order == 2) {
        //    visc = -1*std::pow(2*(expterm(sr/papa_reg) / papa_reg),2)*(p_bg)*inertialNum(sr, p_bg, ro_0, diam, mu_3, A_3, 2*alpha_3);
        //}
    }
    return {visc, include};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Viscosity_VOF(const amrex::Real sr, const amrex::Real dens, const int order, const amrex::Vector<incflo::FLUID_t>& fluid_vof) 
{
    amrex::Real conc = get_concentration(dens,fluid_vof[0].rho,fluid_vof[1].rho);
    auto [visc0, include0] = Viscosity_Single(sr,order,fluid_vof[0]);
    auto [visc1, include1] = Viscosity_Single(sr,order,fluid_vof[1]);

    amrex::Real visc;
    if (include0 and include1) visc = mixture_viscosity(conc,visc0,visc1);
    if (!include0 and !include1) visc = 0.0;
    if (!include0 and include1) visc = visc1;
    if (include0 and !include1) visc = visc0;

    return visc;
}

}

void incflo::compute_viscosity (amrex::Vector<amrex::MultiFab*> const& vel_eta,
                                amrex::Vector<amrex::MultiFab*> const& rho,
                                amrex::Vector<amrex::MultiFab*> const& vel,
                                amrex::Real time, int nghost, int order)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        compute_viscosity_at_level(lev, vel_eta[lev], rho[lev], vel[lev], geom[lev], time, nghost, order);
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
                                         amrex::Geometry& lev_geom,
                                         amrex::Real /*time*/, int nghost, int order)
{

#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

    Real idx = Real(1.0) / lev_geom.CellSize(0);
    Real idy = Real(1.0) / lev_geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
    Real idz = Real(1.0) / lev_geom.CellSize(2);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_eta,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.growntilebox(nghost);
        Array4<Real> const& eta_arr = vel_eta->array(mfi);
        Array4<Real const> const& vel_arr = vel->const_array(mfi);
        Array4<Real const> const& rho_arr = rho->const_array(mfi);
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
                Real sr = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);
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
    }
}

void incflo::compute_tracer_diff_coeff (Vector<MultiFab*> const& tra_eta, int nghost)
{
    for (auto mf : tra_eta) {
        for (int n = 0; n < m_ntrac; ++n) {
            mf->setVal(m_mu_s[n], n, 1, nghost);
        }
    }
}
