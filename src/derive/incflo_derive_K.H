#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate (int i, int j, int k,
                               AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 3)
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    Real wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    Real wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    Real uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    Real vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    Real wz = amrex::Real(0.5) * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                     + (uy+vx)*(uy+vx));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_nodal (int i, int j, int k,
                                     AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                                     amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 3)

    Real uxpp = (vel(i,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real uxpm = (vel(i,j-1,k,0) - vel(i-1,j-1,k,0)) * idx;
    Real uxmp = (vel(i,j,k-1,0) - vel(i-1,j,k-1,0)) * idx;
    Real uxmm = (vel(i,j-1,k-1,0) - vel(i-1,j-1,k-1,0)) * idx;
    Real vxpp = (vel(i,j,k,1) - vel(i-1,j,k,1)) * idx;
    Real vxpm = (vel(i,j-1,k,1) - vel(i-1,j-1,k,1)) * idx;
    Real vxmp = (vel(i,j,k-1,1) - vel(i-1,j,k-1,1)) * idx;
    Real vxmm = (vel(i,j-1,k-1,1) - vel(i-1,j-1,k-1,1)) * idx;
    Real wxpp = (vel(i,j,k,2) - vel(i-1,j,k,2)) * idx;
    Real wxpm = (vel(i,j-1,k,2) - vel(i-1,j-1,k,2)) * idx;
    Real wxmp = (vel(i,j,k-1,2) - vel(i-1,j,k-1,2)) * idx;
    Real wxmm = (vel(i,j-1,k-1,2) - vel(i-1,j-1,k-1,2)) * idx;

    Real uypp = (vel(i,j,k,0) - vel(i,j-1,k,0)) * idy;
    Real uypm = (vel(i-1,j,k,0) - vel(i-1,j-1,k,0)) * idy;
    Real uymp = (vel(i,j,k-1,0) - vel(i,j-1,k-1,0)) * idy;
    Real uymm = (vel(i-1,j,k-1,0) - vel(i-1,j-1,k-1,0)) * idy;
    Real vypp = (vel(i,j,k,1) - vel(i,j-1,k,1)) * idy;
    Real vypm = (vel(i-1,j,k,1) - vel(i-1,j-1,k,1)) * idy;
    Real vymp = (vel(i,j,k-1,1) - vel(i,j-1,k-1,1)) * idy;
    Real vymm = (vel(i-1,j,k-1,1) - vel(i-1,j-1,k-1,1)) * idy;
    Real wypp = (vel(i,j,k,2) - vel(i,j-1,k,2)) * idy;
    Real wypm = (vel(i-1,j,k,2) - vel(i-1,j-1,k,2)) * idy;
    Real wymp = (vel(i,j,k-1,2) - vel(i,j-1,k-1,2)) * idy;
    Real wymm = (vel(i-1,j,k-1,2) - vel(i-1,j-1,k-1,2)) * idy;

    Real uzpp = (vel(i,j,k,0) - vel(i,j,k-1,0)) * idz;
    Real uzpm = (vel(i-1,j,k,0) - vel(i-1,j,k-1,0)) * idz;
    Real uzmp = (vel(i,j-1,k,0) - vel(i,j-1,k-1,0)) * idz;
    Real uzmm = (vel(i-1,j-1,k,0) - vel(i-1,j-1,k-1,0)) * idz;
    Real vzpp = (vel(i,j,k,1) - vel(i,j,k-1,1)) * idz;
    Real vzpm = (vel(i-1,j,k,1) - vel(i-1,j,k-1,1)) * idz;
    Real vzmp = (vel(i,j-1,k,1) - vel(i,j-1,k-1,1)) * idz;
    Real vzmm = (vel(i-1,j-1,k,1) - vel(i-1,j-1,k-1,1)) * idz;
    Real wzpp = (vel(i,j,k,2) - vel(i,j,k-1,2)) * idz;
    Real wzpm = (vel(i-1,j,k,2) - vel(i-1,j,k-1,2)) * idz;
    Real wzmp = (vel(i,j-1,k,2) - vel(i,j-1,k-1,2)) * idz;
    Real wzmm = (vel(i-1,j-1,k,2) - vel(i-1,j-1,k-1,2)) * idz;

    Real ux = amrex::Real(0.25)*(uxpp + uxpm + uxmp + uxmm);
    Real uy = amrex::Real(0.25)*(uypp + uypm + uymp + uymm);
    Real uz = amrex::Real(0.25)*(uzpp + uzpm + uzmp + uzmm);
    Real vx = amrex::Real(0.25)*(vxpp + vxpm + vxmp + vxmm);
    Real vy = amrex::Real(0.25)*(vypp + vypm + vymp + vymm);
    Real vz = amrex::Real(0.25)*(vzpp + vzpm + vzmp + vzmm);
    Real wx = amrex::Real(0.25)*(wxpp + wxpm + wxmp + wxmm);
    Real wy = amrex::Real(0.25)*(wypp + wypm + wymp + wymm);
    Real wz = amrex::Real(0.25)*(wzpp + wzpm + wzmp + wzmm);

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                     + (uy+vx)*(uy+vx));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_rho_nodal (int i, int j, int k, amrex::Array4<amrex::Real const> const& rho) noexcept
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 3)

    return Real(0.125)*(  rho(i-1,j-1,k-1) + rho(i,j-1,k-1)  + 
                          rho(i-1,j,k-1)   + rho(i,j,k-1)    +
                          rho(i-1,j-1,k)   + rho(i,j-1,k)    +
                          rho(i-1,j,k)     + rho(i,j,k) );

#endif
}

#ifdef AMREX_USE_EB
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_eb (int i, int j, int k,
                                  AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                                  amrex::Array4<amrex::Real const> const& vel,
                                  amrex::EBCellFlag flag) noexcept
{
    using namespace amrex;

    if (flag.isRegular()) {
        return incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel);
    } else if (flag.isSingleValued()) {
        constexpr Real c0 = -1.5;
        constexpr Real c1 =  amrex::Real(2.0);
        constexpr Real c2 = -amrex::Real(0.5);
        Real ux, vx, uy, vy;
#if (AMREX_SPACEDIM == 3)
        Real wx, wy, uz, vz, wz;
#endif
        if (!flag.isConnected( 1,0,0))
        {
            // Covered cell to the right, go fish left
            ux = - (c0 * vel(i  ,j,k,0)
                  + c1 * vel(i-1,j,k,0)
                  + c2 * vel(i-2,j,k,0)) * idx;
            vx = - (c0 * vel(i  ,j,k,1)
                  + c1 * vel(i-1,j,k,1)
                  + c2 * vel(i-2,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = - (c0 * vel(i  ,j,k,2)
                  + c1 * vel(i-1,j,k,2)
                  + c2 * vel(i-2,j,k,2)) * idx;
#endif
        }
        else if (!flag.isConnected(-1,0,0))
        {
            // Covered cell to the left, go fish right
            ux = (c0 * vel(i  ,j,k,0)
                + c1 * vel(i+1,j,k,0)
                + c2 * vel(i+2,j,k,0)) * idx;
            vx = (c0 * vel(i  ,j,k,1)
                + c1 * vel(i+1,j,k,1)
                + c2 * vel(i+2,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = (c0 * vel(i  ,j,k,2)
                + c1 * vel(i+1,j,k,2)
                + c2 * vel(i+2,j,k,2)) * idx;
#endif
        }
        else
        {
            // No covered cells right or left, use standard stencil
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
            vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
#if (AMREX_SPACEDIM == 3)
            wx = amrex::Real(0.5) * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;
#endif
        }

        // Do the same in y-direction
        if (!flag.isConnected(0, 1,0))
        {
            uy = - (c0 * vel(i,j  ,k,0)
                  + c1 * vel(i,j-1,k,0)
                  + c2 * vel(i,j-2,k,0)) * idy;
            vy = - (c0 * vel(i,j  ,k,1)
                  + c1 * vel(i,j-1,k,1)
                  + c2 * vel(i,j-2,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = - (c0 * vel(i,j  ,k,2)
                  + c1 * vel(i,j-1,k,2)
                  + c2 * vel(i,j-2,k,2)) * idy;
#endif
        }
        else if (!flag.isConnected(0,-1,0))
        {
            uy = (c0 * vel(i,j  ,k,0)
                + c1 * vel(i,j+1,k,0)
                + c2 * vel(i,j+2,k,0)) * idy;
            vy = (c0 * vel(i,j  ,k,1)
                + c1 * vel(i,j+1,k,1)
                + c2 * vel(i,j+2,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = (c0 * vel(i,j  ,k,2)
                + c1 * vel(i,j+1,k,2)
                + c2 * vel(i,j+2,k,2)) * idy;
#endif
        }
        else
        {
            uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
#if (AMREX_SPACEDIM == 3)
            wy = amrex::Real(0.5) * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;
#endif
        }
#if (AMREX_SPACEDIM == 3)
        // Do the same in z-direction
        if (!flag.isConnected(0,0, 1))
        {
            uz = - (c0 * vel(i,j,k  ,0)
                  + c1 * vel(i,j,k-1,0)
                  + c2 * vel(i,j,k-2,0)) * idz;
            vz = - (c0 * vel(i,j,k  ,1)
                  + c1 * vel(i,j,k-1,1)
                  + c2 * vel(i,j,k-2,1)) * idz;
            wz = - (c0 * vel(i,j,k  ,2)
                  + c1 * vel(i,j,k-1,2)
                  + c2 * vel(i,j,k-2,2)) * idz;
        }
        else if (!flag.isConnected(0,0,-1))
        {
            uz = (c0 * vel(i,j,k  ,0)
                + c1 * vel(i,j,k+1,0)
                + c2 * vel(i,j,k+2,0)) * idz;
            vz = (c0 * vel(i,j,k  ,1)
                + c1 * vel(i,j,k+1,1)
                + c2 * vel(i,j,k+2,1)) * idz;
            wz = (c0 * vel(i,j,k  ,2)
                + c1 * vel(i,j,k+1,2)
                + c2 * vel(i,j,k+2,2)) * idz;
        }
        else
        {
            uz = amrex::Real(0.5) * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
            vz = amrex::Real(0.5) * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
            wz = amrex::Real(0.5) * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;
        }

        return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                         + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
        return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                         + (uy+vx)*(uy+vx));
#endif
    } else {
        return 0.0;
    }
}
#endif

#endif
