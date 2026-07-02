#include <hydro_weno5.H>
#include <hydro_weno5_edge_state_K.H>

using namespace amrex;

void
WENO5::ExtrapVelToFaces(MultiFab const& a_vel,
                        AMREX_D_DECL(MultiFab& a_umac,
                                     MultiFab& a_vmac,
                                     MultiFab& a_wmac),
                        Geometry const& a_geom,
                        Vector<BCRec> const& h_bcrec,
                        BCRec const* d_bcrec,
                        bool allow_inflow_on_outflow)
{
    BL_PROFILE("WENO5::ExtrapVelToFaces");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        AMREX_D_TERM(Box const& ubx = mfi.nodaltilebox(0);,
                     Box const& vbx = mfi.nodaltilebox(1);,
                     Box const& wbx = mfi.nodaltilebox(2););

        AMREX_D_TERM(Array4<Real> const& u = a_umac.array(mfi);,
                     Array4<Real> const& v = a_vmac.array(mfi);,
                     Array4<Real> const& w = a_wmac.array(mfi););

        Array4<Real const> const& vcc = a_vel.const_array(mfi);
        ExtrapVelToFacesBox(AMREX_D_DECL(ubx, vbx, wbx),
                            AMREX_D_DECL(u, v, w),
                            vcc, a_geom, h_bcrec, d_bcrec,
                            allow_inflow_on_outflow);
    }
}

void
WENO5::ExtrapVelToFacesBox(AMREX_D_DECL(Box const& ubx,
                                        Box const& vbx,
                                        Box const& wbx),
                           AMREX_D_DECL(Array4<Real> const& u,
                                        Array4<Real> const& v,
                                        Array4<Real> const& w),
                           Array4<Real const> const& vcc,
                           Geometry const& geom,
                           Vector<BCRec> const& /*h_bcrec*/,
                           BCRec const* d_bcrec,
                           bool allow_inflow_on_outflow)
{
    BL_PROFILE("WENO5::ExtrapVelToFacesBox");

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        u(i,j,k) = WENO5::hydro_weno5_xface_velocity(i, j, k, vcc, d_bcrec,
                                                     dlo.x, dhi.x,
                                                     allow_inflow_on_outflow);
    });

    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        v(i,j,k) = WENO5::hydro_weno5_yface_velocity(i, j, k, vcc, d_bcrec,
                                                     dlo.y, dhi.y,
                                                     allow_inflow_on_outflow);
    });

#if (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(wbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        w(i,j,k) = WENO5::hydro_weno5_zface_velocity(i, j, k, vcc, d_bcrec,
                                                     dlo.z, dhi.z,
                                                     allow_inflow_on_outflow);
    });
#endif
}
