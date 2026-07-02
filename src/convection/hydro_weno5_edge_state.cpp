#include <hydro_weno5.H>
#include <hydro_weno5_edge_state_K.H>
#include <hydro_bcs_K.H>

using namespace amrex;

void
WENO5::ComputeEdgeState(Box const& bx,
                        AMREX_D_DECL(Array4<Real> const& xedge,
                                     Array4<Real> const& yedge,
                                     Array4<Real> const& zedge),
                        Array4<Real const> const& q,
                        int ncomp,
                        AMREX_D_DECL(Array4<Real const> const& umac,
                                     Array4<Real const> const& vmac,
                                     Array4<Real const> const& wmac),
                        Box const& domain,
                        Vector<BCRec> const& /*bcs*/,
                        BCRec const* d_bcrec,
                        bool is_velocity,
                        Array4<int const> const& bc_arr)
{
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    AMREX_D_TERM(Box const& xbx = amrex::surroundingNodes(bx, 0);,
                 Box const& ybx = amrex::surroundingNodes(bx, 1);,
                 Box const& zbx = amrex::surroundingNodes(bx, 2););

    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, d_bcrec, bc_arr);
        xedge(i,j,k,n) = WENO5::hydro_weno5_xedge_state(i, j, k, n, q, umac,
                                                        bc, dlo.x, dhi.x,
                                                        is_velocity);
    });

    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, d_bcrec, bc_arr);
        yedge(i,j,k,n) = WENO5::hydro_weno5_yedge_state(i, j, k, n, q, vmac,
                                                        bc, dlo.y, dhi.y,
                                                        is_velocity);
    });

#if (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(zbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = HydroBC::getBC(i, j, k, n, domain, d_bcrec, bc_arr);
        zedge(i,j,k,n) = WENO5::hydro_weno5_zedge_state(i, j, k, n, q, wmac,
                                                        bc, dlo.z, dhi.z,
                                                        is_velocity);
    });
#endif
}
