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
    amrex::Real visc;
    if (conc < 1.0e-12) visc = visc1;
    else if (conc > amrex::Real(1.0) - 1.0e-12) visc = visc0;
    else visc = 1.0/((conc/visc0) + ((1.0-conc)/visc1));
    return visc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::tuple<amrex::Real, bool> Viscosity_Single(const amrex::Real sr, const int order, const incflo::FLUID_t& fluid, const amrex::Real press, const bool slip_surface) 
{
    amrex::Real visc = 0.0;
    bool include = false;
    
    if (fluid.fluid_model == incflo::FluidModel::Newtonian) {
        if (order == 0) {visc = fluid.mu; include = true;}
        else {visc = 0.0; include = false;}
    }
    else if (fluid.fluid_model == incflo::FluidModel::Powerlaw) {
        if (order == 0) {
            if (sr > Real(1.0e-13)) {
                visc = (fluid.mu * std::pow(sr,fluid.n_0-1.0));
            }
            else {
                visc = fluid.mu;
            }
            include = true;
        }
        else if (order == 1) {
            if (sr > Real(1.0e-13)) {
                visc = (fluid.mu_1 * std::pow(sr,fluid.n_1-1.0));
            }
            else {
                visc = fluid.mu;
            }
            include = true;
        }
    }
    else if (fluid.fluid_model == incflo::FluidModel::Bingham) {
        
        if (order == 0) {
            amrex::Real compute_visc;
            if ((fluid.papa_reg > Real(0.0)) and (fluid.max_visc > Real(0.0))) amrex::Abort("can not prescribe both regularization parameter and maximum viscosity");
            
            if (fluid.papa_reg > Real(0.0)) { // Papanastasiou regularization
                visc = fluid.mu + fluid.tau_0 * expterm(sr/fluid.papa_reg) / fluid.papa_reg;
            }
            else if (fluid.max_visc > Real(0.0)) { // ceiling viscosity regularization
                if (sr > Real(1.0e-13)) {
                    compute_visc = fluid.mu + fluid.tau_0/sr;
                    visc = std::min(compute_visc, fluid.max_visc);
                }
                else {
                    visc = fluid.max_visc;
                }
            }
            include = true;
        }
        else if (order == 1) {
            amrex::Real compute_visc;
            if ((fluid.papa_reg_1 > Real(0.0)) and (fluid.max_visc_1 > Real(0.0))) amrex::Abort("can not prescribe both regularization parameter and maximum viscosity");
            
            if (fluid.papa_reg_1 > Real(0.0)) { // Papanastasiou regularization
                visc = fluid.mu_1 + fluid.tau_1 * expterm(sr/fluid.papa_reg_1) / fluid.papa_reg_1;
            }
            else if (fluid.max_visc_1 > Real(0.0)) { // ceiling viscosity regularization
                if (sr > Real(1.0e-13)) {
                    compute_visc = fluid.mu_1 + fluid.tau_1/sr;
                    visc = std::min(compute_visc, fluid.max_visc_1);
                }
                else {
                    visc = fluid.max_visc_1;
                }
            }
            include = true;
        }
        return {visc, include};

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
            
            amrex::Real compute_visc;

            if ((fluid.papa_reg > Real(0.0)) and (fluid.max_visc > Real(0.0))) amrex::Abort("can not prescribe both regularization parameter and maximum viscosity");

            if (press < 0.0) { // negative pressure
                visc = fluid.min_visc;
                include = true;
            }
            else { // positive pressure
                if (fluid.papa_reg > Real(0.0)) { // Papanastasiou regularization
                    compute_visc = (fluid.tau_0*press)*expterm(sr/fluid.papa_reg)/fluid.papa_reg + 
                                    fluid.A_0*std::pow(fluid.diam,fluid.alpha_0)*std::pow(fluid.rho,0.5*fluid.alpha_0)*
                                    std::pow(press,1.0-0.5*fluid.alpha_0)*std::pow(sr,fluid.alpha_0); // regularized viscosity
                }
                else if (fluid.max_visc > Real(0.0)) { // ceiling viscosity regularization
                    if (sr > Real(1.0e-13)) {
                        amrex::Real compute_visc_org = (fluid.tau_0*press + 
                                                        fluid.A_0*std::pow(fluid.diam,fluid.alpha_0)*std::pow(fluid.rho,0.5*fluid.alpha_0)*
                                                        std::pow(press,1.0-0.5*fluid.alpha_0)*std::pow(sr,fluid.alpha_0))/sr; // unregularized viscosity
                        compute_visc = std::min(compute_visc_org, fluid.max_visc);
                    }
                    else {
                        compute_visc = fluid.max_visc;
                    }

                    // we will set viscosity to its minimum value if it is a slip-surface (\dot{\gamma} = 0)
                    if (slip_surface) compute_visc = fluid.min_visc;
                }
                visc = std::max(fluid.min_visc, compute_visc); // floor viscosity
                include = true;

            }

        }
        else if (order == 1) {

            amrex::Real compute_visc;

//            if ((fluid.max_visc_1 < Real(0.0)) or (fluid.min_visc_1 < Real(0.0)))
//                amrex::Abort("need positive max and min firs notmal viscosity");

            if (press < 0.0) { // negative pressure
                visc = fluid.min_visc_1;
                include = true;
            }
            else { // positive pressure
                if (sr > Real(1.0e-13)) {
                    amrex::Real compute_visc_org =
                        (1.0/sr/sr)*(fluid.tau_1*press +
                        fluid.A_1*std::pow(fluid.diam,2.0*fluid.alpha_1)*std::pow(fluid.rho,fluid.alpha_1)*
                        std::pow(press,1.0-fluid.alpha_1)*std::pow(sr,2.0*fluid.alpha_1)); // unregularized viscosity

                    compute_visc = std::min(compute_visc_org, fluid.max_visc_1);
                }
                else {
                    compute_visc = fluid.max_visc_1;
                }
                visc = std::max(fluid.min_visc_1, compute_visc);
                include = true;
            }
        }
        //else if (order == 2) {
        //    visc = -1*std::pow(2*(expterm(sr/papa_reg) / papa_reg),2)*(p_bg)*inertialNum(sr, p_bg, ro_0, diam, mu_3, A_3, 2*alpha_3);
        //}
    }
    return {visc, include};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Viscosity_VOF(const amrex::Real sr, const amrex::Real conc, const int order, const amrex::Vector<incflo::FLUID_t>& fluid_vof, const amrex::Real press, bool slip_surface) 
{
//    amrex::Real conc = get_concentration(dens,fluid_vof[0].rho,fluid_vof[1].rho);
    auto [visc0, include0] = Viscosity_Single(sr,order,fluid_vof[0],press,slip_surface);
    auto [visc1, include1] = Viscosity_Single(sr,order,fluid_vof[1],press,slip_surface);

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

bool is_slip_surface(const int i, const int j, const int k, 
                     const amrex::Dim3& domlo, const amrex::Dim3& domhi,
                     amrex::Vector<amrex::Array<int,2>> bc_type)
{
    bool slip_surface = false;
    if ((i<=domlo.x) && (bc_type[0][0] == 1)) slip_surface = true;
    if ((i>=domhi.x) && (bc_type[0][1] == 1)) slip_surface = true;
    if ((j<=domlo.y) && (bc_type[1][0] == 1)) slip_surface = true;
    if ((j>=domhi.y) && (bc_type[1][1] == 1)) slip_surface = true;
#if (AMREX_SPACEDIM==3)
    if ((k<=domlo.z) && (bc_type[2][0] == 1)) slip_surface = true;
    if ((k>=domhi.z) && (bc_type[2][1] == 1)) slip_surface = true;
#endif
    return slip_surface;
}

}

void incflo::compute_viscosity (amrex::Vector<amrex::MultiFab      *> const& vel_eta,
                                amrex::Vector<amrex::MultiFab const*> const& rho,
                                amrex::Vector<amrex::MultiFab const*> const& vel,
                                amrex::Vector<amrex::MultiFab const*> const& press,
                                amrex::Vector<amrex::MultiFab const*> const& press0,
                                amrex::Vector<amrex::MultiFab      *> const& press_visc,
                                amrex::Real time, int nghost, int order)
{
    

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        compute_viscosity_at_level(lev, vel_eta[lev], rho[lev], vel[lev], press[lev], press0[lev], press_visc[lev], geom[lev], time, nghost, order);
    }
}

#ifdef AMREX_USE_EB
void incflo::compute_viscosity_at_level (int lev,
#else
void incflo::compute_viscosity_at_level (int /*lev*/,
#endif
                                         amrex::MultiFab      * vel_eta,
                                         amrex::MultiFab const* rho,
                                         amrex::MultiFab const* vel,
                                         amrex::MultiFab const* press,
                                         amrex::MultiFab const* press0,
                                         amrex::MultiFab      * press_visc,
                                         amrex::Geometry& lev_geom,
                                         amrex::Real /*time*/, int /*nghost*/, int order)
{

    // Get pressure at the top of the domain
    amrex::Gpu::DeviceVector<Real> press_surface;
#if (AMREX_SPACEDIM == 3)
        press_surface.resize((m_n_cells[0]+1)*(m_n_cells[1]+1), -1.0e99);
        getsurfacepressure(press,press_surface,lev_geom);
        ParallelDescriptor::ReduceRealMax(press_surface.data(),(m_n_cells[0]+1)*(m_n_cells[1]+1));
#elif (AMREX_SPACEDIM == 2)
        press_surface.resize(m_n_cells[0]+1, -1.0e99);
        getsurfacepressure(press,press_surface,lev_geom);
        ParallelDescriptor::ReduceRealMax(press_surface.data(),m_n_cells[0]+1);
#endif

#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif
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

        Array4<Real      > const& eta_arr       = vel_eta->array(mfi);
        Array4<Real const> const& vel_arr       = vel->const_array(mfi);
        Array4<Real const> const& rho_arr       = rho->const_array(mfi);
        Array4<Real const> const& p_arr         = press->const_array(mfi);
        Array4<Real const> const& p0_arr        = press0->const_array(mfi);
        Array4<Real      > const& p_visc_arr    = press_visc->array(mfi);
        
        Real* p_surface = press_surface.data();

#ifdef AMREX_USE_EB
        auto const& flag_fab = flags[mfi];
        auto typ = flag_fab.getType(bx);
        if (typ == FabType::covered)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                eta_arr(i,j,k,0) = Real(0.0);
                eta_arr(i,j,k,1) = Real(0.0);
                eta_arr(i,j,k,2) = Real(0.0);
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

                // Get Strain Rate (nodal)
                Real sr = 0.0;
                if (m_need_strainrate) {
                    sr = incflo_strainrate_nodal(i,j,k,AMREX_D_DECL(idx,idy,idz),
                                                  vel_arr,domlo,domhi,bc_type);
                }

                // Get concentration (nodal)
                Real conc;
//                Real dens = m_ro_0;                
                if (m_do_vof) {
                    conc = incflo_nodal_conc(i,j,k,rho_arr,m_fluid_vof[0].rho,m_fluid_vof[1].rho);
                    //dens = incflo_cc_to_nodal(i,j,k,rho_arr,true);
                }

                // Get Pressure (nodal)
                Real pressure = 0.0;
                if (m_p_dep_visc == 1) {
#if (AMREX_SPACEDIM == 3)
                    int index = j*(m_n_cells[0]+1) + i;
                    Real p_surface_ij = p_surface[index];
                    pressure = p_arr(i,j,k) - p_surface_ij + m_p_top_surface;
#elif (AMREX_SPACEDIM == 2)
                    int index = i;
                    Real p_surface_i = p_surface[index];
                    pressure = p_arr(i,j,k) - p_surface_i + m_p_top_surface;
#endif
                }
                else if (m_p_dep_visc == 0) {
                    pressure = p0_arr(i,j,k);
                }
                else {
                    pressure = p_arr(i,j,k);
                }
                p_visc_arr(i,j,k) = pressure;
                
                // Is slip-surface? (need for granular)
                bool slip_surface = is_slip_surface(i,j,k,domlo,domhi,bc_type);

                // Compute viscosity (nodal)
                if (m_do_vof) {
                    eta_arr(i,j,k) = Viscosity_VOF(sr,conc,order,m_fluid_vof,pressure,slip_surface);
                }
                else {
                    auto [visc, include] = Viscosity_Single(sr,order,m_fluid,pressure,slip_surface);
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


void incflo::getsurfacepressure(amrex::MultiFab const* press, amrex::Gpu::DeviceVector<Real>& p_surface_in, amrex::Geometry& geom)
{
    auto const& problo = geom.ProbLoArray();
    auto const& probhi = geom.ProbHiArray();
    AMREX_D_TERM(const Real l_dx = geom.CellSize(0);,
                 const Real l_dy = geom.CellSize(1);,
                 const Real l_dz = geom.CellSize(2););

    for (MFIter mfi(*press,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real const> const& p_arr = press->const_array(mfi);
        Real* p_surface = p_surface_in.data();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 3)
            Real z = problo[2] + k*l_dz;
            if ((z<probhi[2]+Real(0.01)*l_dz) && (z>probhi[2]-Real(0.01)*l_dz)) { // top surface
                int index = j*(m_n_cells[0]+1) + i;
                p_surface[index] = p_arr(i,j,k);
            }
#elif (AMREX_SPACEDIM == 2)
            Real y = problo[1] + j*l_dy;
            if ((y<probhi[1]+Real(0.01)*l_dy) && (y>probhi[1]-Real(0.01)*l_dy)) { // top surface
                int index = i;
                p_surface[index] = p_arr(i,j,k);
            }
#endif
        }); // end MFITer 
    }
}
