#ifndef HYDRO_WENO5_EDGE_STATE_K_H_
#define HYDRO_WENO5_EDGE_STATE_K_H_

#include <hydro_bcs_K.H>
#include <hydro_constants.H>
#include <hydro_godunov_ppm.H>

namespace WENO5 {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_select(amrex::Real qmns, amrex::Real qpls, amrex::Real vel) noexcept
{
    if (vel > small_vel) {
        return qmns;
    } else if (vel < -small_vel) {
        return qpls;
    } else {
        return amrex::Real(0.5) * (qmns + qpls);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_xmns(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge2(q(i-3,j,k,n), q(i-2,j,k,n), q(i-1,j,k,n),
                                q(i  ,j,k,n), q(i+1,j,k,n));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_xpls(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge1(q(i-2,j,k,n), q(i-1,j,k,n), q(i  ,j,k,n),
                                q(i+1,j,k,n), q(i+2,j,k,n));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_ymns(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge2(q(i,j-3,k,n), q(i,j-2,k,n), q(i,j-1,k,n),
                                q(i,j  ,k,n), q(i,j+1,k,n));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_ypls(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge1(q(i,j-2,k,n), q(i,j-1,k,n), q(i,j  ,k,n),
                                q(i,j+1,k,n), q(i,j+2,k,n));
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_zmns(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge2(q(i,j,k-3,n), q(i,j,k-2,n), q(i,j,k-1,n),
                                q(i,j,k  ,n), q(i,j,k+1,n));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5_zpls(int i, int j, int k, int n,
                       amrex::Array4<amrex::Real const> const& q) noexcept
{
    return PPM::weno_js::sedge1(q(i,j,k-2,n), q(i,j,k-1,n), q(i,j,k  ,n),
                                q(i,j,k+1,n), q(i,j,k+2,n));
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_xedge_state(int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& umac,
                                    amrex::BCRec const& bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    if (bcrec.lo(0) == amrex::BCType::ext_dir && i <= domlo) {
        return q(domlo-1,j,k,n);
    }
    if (bcrec.hi(0) == amrex::BCType::ext_dir && i >= domhi+1) {
        return q(domhi+1,j,k,n);
    }

    amrex::Real qmns = weno5_xmns(i, j, k, n, q);
    amrex::Real qpls = weno5_xpls(i, j, k, n, q);

    HydroBC::SetEdgeBCsLo(0, i, j, k, n, q, qmns, qpls, umac(i,j,k),
                          bcrec.lo(0), domlo, is_velocity);
    HydroBC::SetEdgeBCsHi(0, i, j, k, n, q, qmns, qpls, umac(i,j,k),
                          bcrec.hi(0), domhi, is_velocity);

    if ((i == domlo) && (bcrec.lo(0) == amrex::BCType::foextrap ||
                         bcrec.lo(0) == amrex::BCType::hoextrap)) {
        if (umac(i,j,k) >= amrex::Real(0.0) && n == XVEL && is_velocity) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
        }
        qmns = qpls;
    }
    if ((i == domhi+1) && (bcrec.hi(0) == amrex::BCType::foextrap ||
                           bcrec.hi(0) == amrex::BCType::hoextrap)) {
        if (umac(i,j,k) <= amrex::Real(0.0) && n == XVEL && is_velocity) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
        }
        qpls = qmns;
    }

    return weno5_select(qmns, qpls, umac(i,j,k)); // this is right (also for y, z), but just make sure the epsilon is small
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_yedge_state(int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::BCRec const& bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    if (bcrec.lo(1) == amrex::BCType::ext_dir && j <= domlo) {
        return q(i,domlo-1,k,n);
    }
    if (bcrec.hi(1) == amrex::BCType::ext_dir && j >= domhi+1) {
        return q(i,domhi+1,k,n);
    }

    amrex::Real qmns = weno5_ymns(i, j, k, n, q);
    amrex::Real qpls = weno5_ypls(i, j, k, n, q);

    HydroBC::SetEdgeBCsLo(1, i, j, k, n, q, qmns, qpls, vmac(i,j,k),
                          bcrec.lo(1), domlo, is_velocity);
    HydroBC::SetEdgeBCsHi(1, i, j, k, n, q, qmns, qpls, vmac(i,j,k),
                          bcrec.hi(1), domhi, is_velocity);

    if ((j == domlo) && (bcrec.lo(1) == amrex::BCType::foextrap ||
                         bcrec.lo(1) == amrex::BCType::hoextrap)) {
        if (vmac(i,j,k) >= amrex::Real(0.0) && n == YVEL && is_velocity) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
        }
        qmns = qpls;
    }
    if ((j == domhi+1) && (bcrec.hi(1) == amrex::BCType::foextrap ||
                           bcrec.hi(1) == amrex::BCType::hoextrap)) {
        if (vmac(i,j,k) <= amrex::Real(0.0) && n == YVEL && is_velocity) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
        }
        qpls = qmns;
    }

    return weno5_select(qmns, qpls, vmac(i,j,k));
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_zedge_state(int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& wmac,
                                    amrex::BCRec const& bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    if (bcrec.lo(2) == amrex::BCType::ext_dir && k <= domlo) {
        return q(i,j,domlo-1,n);
    }
    if (bcrec.hi(2) == amrex::BCType::ext_dir && k >= domhi+1) {
        return q(i,j,domhi+1,n);
    }

    amrex::Real qmns = weno5_zmns(i, j, k, n, q);
    amrex::Real qpls = weno5_zpls(i, j, k, n, q);

    HydroBC::SetEdgeBCsLo(2, i, j, k, n, q, qmns, qpls, wmac(i,j,k),
                          bcrec.lo(2), domlo, is_velocity);
    HydroBC::SetEdgeBCsHi(2, i, j, k, n, q, qmns, qpls, wmac(i,j,k),
                          bcrec.hi(2), domhi, is_velocity);

    if ((k == domlo) && (bcrec.lo(2) == amrex::BCType::foextrap ||
                         bcrec.lo(2) == amrex::BCType::hoextrap)) {
        if (wmac(i,j,k) >= amrex::Real(0.0) && n == ZVEL && is_velocity) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
        }
        qmns = qpls;
    }
    if ((k == domhi+1) && (bcrec.hi(2) == amrex::BCType::foextrap ||
                           bcrec.hi(2) == amrex::BCType::hoextrap)) {
        if (wmac(i,j,k) <= amrex::Real(0.0) && n == ZVEL && is_velocity) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
        }
        qpls = qmns;
    }

    return weno5_select(qmns, qpls, wmac(i,j,k));
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_xface_velocity(int i, int j, int k,
                                       amrex::Array4<amrex::Real const> const& vcc,
                                       amrex::BCRec const* d_bcrec,
                                       int domlo, int domhi,
                                       bool allow_inflow_on_outflow) noexcept
{
    amrex::Real qmns = weno5_xmns(i, j, k, 0, vcc);
    amrex::Real qpls = weno5_xpls(i, j, k, 0, vcc);

    HydroBC::SetExtrapVelBCsLo(0, i, j, k, 0, vcc, qmns, qpls, d_bcrec[0].lo(0), domlo);
    HydroBC::SetExtrapVelBCsHi(0, i, j, k, 0, vcc, qmns, qpls, d_bcrec[0].hi(0), domhi);

    if (!allow_inflow_on_outflow) {
        if ((i == domlo) && (d_bcrec[0].lo(0) == amrex::BCType::foextrap ||
                             d_bcrec[0].lo(0) == amrex::BCType::hoextrap)) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
            qmns = qpls;
        }
        if ((i == domhi+1) && (d_bcrec[0].hi(0) == amrex::BCType::foextrap ||
                               d_bcrec[0].hi(0) == amrex::BCType::hoextrap)) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
            qpls = qmns;
        }
    }

    if (i == domlo && d_bcrec[0].lo(0) == amrex::BCType::ext_dir) {
        return vcc(i-1,j,k,0);
    }
    if (i == domhi+1 && d_bcrec[0].hi(0) == amrex::BCType::ext_dir) {
        return vcc(i,j,k,0);
    }
    return weno5_select(qmns, qpls, amrex::Real(0.5) * (qmns + qpls)); // this needs to be correct upwinding logic -- Riemann problem (also for y, z)
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_yface_velocity(int i, int j, int k,
                                       amrex::Array4<amrex::Real const> const& vcc,
                                       amrex::BCRec const* d_bcrec,
                                       int domlo, int domhi,
                                       bool allow_inflow_on_outflow) noexcept
{
    amrex::Real qmns = weno5_ymns(i, j, k, 1, vcc);
    amrex::Real qpls = weno5_ypls(i, j, k, 1, vcc);

    HydroBC::SetExtrapVelBCsLo(1, i, j, k, 1, vcc, qmns, qpls, d_bcrec[1].lo(1), domlo);
    HydroBC::SetExtrapVelBCsHi(1, i, j, k, 1, vcc, qmns, qpls, d_bcrec[1].hi(1), domhi);

    if (!allow_inflow_on_outflow) {
        if ((j == domlo) && (d_bcrec[1].lo(1) == amrex::BCType::foextrap ||
                             d_bcrec[1].lo(1) == amrex::BCType::hoextrap)) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
            qmns = qpls;
        }
        if ((j == domhi+1) && (d_bcrec[1].hi(1) == amrex::BCType::foextrap ||
                               d_bcrec[1].hi(1) == amrex::BCType::hoextrap)) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
            qpls = qmns;
        }
    }

    if (j == domlo && d_bcrec[1].lo(1) == amrex::BCType::ext_dir) {
        return vcc(i,j-1,k,1);
    }
    if (j == domhi+1 && d_bcrec[1].hi(1) == amrex::BCType::ext_dir) {
        return vcc(i,j,k,1);
    }
    return weno5_select(qmns, qpls, amrex::Real(0.5) * (qmns + qpls));
}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_weno5_zface_velocity(int i, int j, int k,
                                       amrex::Array4<amrex::Real const> const& vcc,
                                       amrex::BCRec const* d_bcrec,
                                       int domlo, int domhi,
                                       bool allow_inflow_on_outflow) noexcept
{
    amrex::Real qmns = weno5_zmns(i, j, k, 2, vcc);
    amrex::Real qpls = weno5_zpls(i, j, k, 2, vcc);

    HydroBC::SetExtrapVelBCsLo(2, i, j, k, 2, vcc, qmns, qpls, d_bcrec[2].lo(2), domlo);
    HydroBC::SetExtrapVelBCsHi(2, i, j, k, 2, vcc, qmns, qpls, d_bcrec[2].hi(2), domhi);

    if (!allow_inflow_on_outflow) {
        if ((k == domlo) && (d_bcrec[2].lo(2) == amrex::BCType::foextrap ||
                             d_bcrec[2].lo(2) == amrex::BCType::hoextrap)) {
            qpls = amrex::min(qpls, amrex::Real(0.0));
            qmns = qpls;
        }
        if ((k == domhi+1) && (d_bcrec[2].hi(2) == amrex::BCType::foextrap ||
                               d_bcrec[2].hi(2) == amrex::BCType::hoextrap)) {
            qmns = amrex::max(qmns, amrex::Real(0.0));
            qpls = qmns;
        }
    }

    if (k == domlo && d_bcrec[2].lo(2) == amrex::BCType::ext_dir) {
        return vcc(i,j,k-1,2);
    }
    if (k == domhi+1 && d_bcrec[2].hi(2) == amrex::BCType::ext_dir) {
        return vcc(i,j,k,2);
    }
    return weno5_select(qmns, qpls, amrex::Real(0.5) * (qmns + qpls));
}
#endif

}

#endif
