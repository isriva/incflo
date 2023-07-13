#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate (int i, int j, int k,
                               AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                               amrex::Array4<amrex::Real const> const& vel, const amrex::Dim3& domlo, const amrex::Dim3& domhi,
                               amrex::Vector<amrex::Array<int,2>> bc_type) noexcept
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

    // x-wall (lo)
    if (i==domlo.x) {
        if (bc_type[0][0] == 2) { // no-slip
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) + vel(i,j,k,0)) * idx;
            vx = amrex::Real(0.5) * (vel(i+1,j,k,1) + vel(i,j,k,1)) * idx;
            wx = amrex::Real(0.5) * (vel(i+1,j,k,2) + vel(i,j,k,2)) * idx;
        }
        else if (bc_type[0][0] == 1) { // slip
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) + vel(i,j,k,0)) * idx;
        }
    }

    // x-wall (hi)
    if (i==domhi.x) {
        if (bc_type[0][1] == 2) { // no-slip
            ux = amrex::Real(-0.5) * (vel(i-1,j,k,0) + vel(i,j,k,0)) * idx;
            vx = amrex::Real(-0.5) * (vel(i-1,j,k,1) + vel(i,j,k,1)) * idx;
            wx = amrex::Real(-0.5) * (vel(i-1,j,k,2) + vel(i,j,k,2)) * idx;
        }
        else if (bc_type[0][1] == 1) { // slip
            ux = amrex::Real(-0.5) * (vel(i-1,j,k,0) + vel(i,j,k,0)) * idx;
        }
    }

    // y-wall (lo)
    if (j==domlo.y) {
        if (bc_type[1][0] == 2) { // no-slip
            uy = amrex::Real(0.5) * (vel(i,j+1,k,0) + vel(i,j,k,0)) * idy;
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) + vel(i,j,k,1)) * idy;
            wy = amrex::Real(0.5) * (vel(i,j+1,k,2) + vel(i,j,k,2)) * idy;
        }
        else if (bc_type[1][0] == 1) { // slip
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) + vel(i,j,k,1)) * idy;
        }
    }

    // y-wall (hi)
    if (j==domhi.y) {
        if (bc_type[1][1] == 2) { // no-slip
            uy = amrex::Real(-0.5) * (vel(i,j-1,k,0) + vel(i,j,k,0)) * idy;
            vy = amrex::Real(-0.5) * (vel(i,j-1,k,1) + vel(i,j,k,1)) * idy;
            wy = amrex::Real(-0.5) * (vel(i,j-1,k,2) + vel(i,j,k,2)) * idy;
        }
        else if (bc_type[1][1] == 1) { // slip
            vy = amrex::Real(-0.5) * (vel(i,j-1,k,1) + vel(i,j,k,1)) * idy;
        }
    }

    // z-wall (lo)
    if (k==domlo.z) {
        if (bc_type[2][0] == 2) { // no-slip
            uz = amrex::Real(0.5) * (vel(i,j,k+1,0) + vel(i,j,k,0)) * idz;
            vz = amrex::Real(0.5) * (vel(i,j,k+1,1) + vel(i,j,k,1)) * idz;
            wz = amrex::Real(0.5) * (vel(i,j,k+1,2) + vel(i,j,k,2)) * idz;
        }
        else if (bc_type[2][0] == 1) { // slip
            wz = amrex::Real(0.5) * (vel(i,j,k+1,2) + vel(i,j,k,2)) * idz;
        }
    }

    // z-wall (hi)
    if (k==domhi.z) {
        if (bc_type[2][1] == 2) { // no-slip
            uz = amrex::Real(-0.5) * (vel(i,j,k-1,0) + vel(i,j,k,0)) * idz;
            vz = amrex::Real(-0.5) * (vel(i,j,k-1,1) + vel(i,j,k,1)) * idz;
            wz = amrex::Real(-0.5) * (vel(i,j,k-1,2) + vel(i,j,k,2)) * idz;
        }
        else if (bc_type[2][1] == 1) { // slip
            wz = amrex::Real(-0.5) * (vel(i,j,k-1,2) + vel(i,j,k,2)) * idz;
        }
    }

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy + amrex::Real(2.0) * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
#else
    Real ux = amrex::Real(0.5) * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real vx = amrex::Real(0.5) * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;

    Real uy = amrex::Real(0.5) * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    Real vy = amrex::Real(0.5) * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;

    // x-wall (lo)
    if (i==domlo.x) {
        if (bc_type[0][0] == 2) { // no-slip
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) + vel(i,j,k,0)) * idx;
            vx = amrex::Real(0.5) * (vel(i+1,j,k,1) + vel(i,j,k,1)) * idx;
        }
        else if (bc_type[0][0] == 1) { // slip
            ux = amrex::Real(0.5) * (vel(i+1,j,k,0) + vel(i,j,k,0)) * idx;
        }
    }

    // x-wall (hi)
    if (i==domhi.x) {
        if (bc_type[0][1] == 2) { // no-slip
            ux = amrex::Real(-0.5) * (vel(i-1,j,k,0) + vel(i,j,k,0)) * idx;
            vx = amrex::Real(-0.5) * (vel(i-1,j,k,1) + vel(i,j,k,1)) * idx;
        }
        else if (bc_type[0][1] == 1) { // slip
            ux = amrex::Real(-0.5) * (vel(i-1,j,k,0) + vel(i,j,k,0)) * idx;
        }
    }

    // y-wall (lo)
    if (j==domlo.y) {
        if (bc_type[1][0] == 2) { // no-slip
            uy = amrex::Real(0.5) * (vel(i,j+1,k,0) + vel(i,j,k,0)) * idy;
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) + vel(i,j,k,1)) * idy;
        }
        else if (bc_type[1][0] == 1) { // slip
            vy = amrex::Real(0.5) * (vel(i,j+1,k,1) + vel(i,j,k,1)) * idy;
        }
    }

    // y-wall (hi)
    if (j==domhi.y) {
        if (bc_type[1][1] == 2) { // no-slip
            uy = amrex::Real(-0.5) * (vel(i,j-1,k,0) + vel(i,j,k,0)) * idy;
            vy = amrex::Real(-0.5) * (vel(i,j-1,k,1) + vel(i,j,k,1)) * idy;
        }
        else if (bc_type[1][1] == 1) { // slip
            vy = amrex::Real(-0.5) * (vel(i,j-1,k,1) + vel(i,j,k,1)) * idy;
        }
    }

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                     + (uy+vx)*(uy+vx));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_nodal (int i, int j, int k,
                                     AMREX_D_DECL(amrex::Real idx, amrex::Real idy, amrex::Real idz),
                                     amrex::Array4<amrex::Real const> const& vel, const amrex::Dim3& domlo, const amrex::Dim3& domhi,
                                     amrex::Vector<amrex::Array<int,2>> bc_type) noexcept
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

    // x-wall (lo)
    if (i==domlo.x) {
        if (bc_type[0][0] == 2) { // no-slip
            uxpp = amrex::Real(2.0) * (vel(i,j,k,0)) * idx;
            uxpm = amrex::Real(2.0) * (vel(i,j-1,k,0)) * idx;
            uxmp = amrex::Real(2.0) * (vel(i,j,k-1,0)) * idx;
            uxmm = amrex::Real(2.0) * (vel(i,j-1,k-1,0)) * idx;
            vxpp = amrex::Real(2.0) * (vel(i,j,k,1)) * idx;
            vxpm = amrex::Real(2.0) * (vel(i,j-1,k,1)) * idx;
            vxmp = amrex::Real(2.0) * (vel(i,j,k-1,1)) * idx;
            vxmm = amrex::Real(2.0) * (vel(i,j-1,k-1,1)) * idx;
            wxpp = amrex::Real(2.0) * (vel(i,j,k,2)) * idx;
            wxpm = amrex::Real(2.0) * (vel(i,j-1,k,2)) * idx;
            wxmp = amrex::Real(2.0) * (vel(i,j,k-1,2)) * idx;
            wxmm = amrex::Real(2.0) * (vel(i,j-1,k-1,2)) * idx;
        }
        else if (bc_type[0][0] == 1) { // slip
            uxpp = amrex::Real(2.0) * (vel(i,j,k,0)) * idx;
            uxpm = amrex::Real(2.0) * (vel(i,j-1,k,0)) * idx;
            uxmp = amrex::Real(2.0) * (vel(i,j,k-1,0)) * idx;
            uxmm = amrex::Real(2.0) * (vel(i,j-1,k-1,0)) * idx;
            vxpp = amrex::Real(0.0);
            vxpm = amrex::Real(0.0);
            vxmp = amrex::Real(0.0);
            vxmm = amrex::Real(0.0);
            wxpp = amrex::Real(0.0);
            wxpm = amrex::Real(0.0);
            wxmp = amrex::Real(0.0);
            wxmm = amrex::Real(0.0);
        }
    }
    // x-wall (hi)
    if (i==domhi.x + 1) {
        if (bc_type[0][1] == 2) { // no-slip
            uxpp = amrex::Real(-2.0) * (vel(i-1,j,k,0)) * idx;
            uxpm = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idx;
            uxmp = amrex::Real(-2.0) * (vel(i-1,j,k-1,0)) * idx;
            uxmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,0)) * idx;
            vxpp = amrex::Real(-2.0) * (vel(i-1,j,k,1)) * idx;
            vxpm = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idx;
            vxmp = amrex::Real(-2.0) * (vel(i-1,j,k-1,1)) * idx;
            vxmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,1)) * idx;
            wxpp = amrex::Real(-2.0) * (vel(i-1,j,k,2)) * idx;
            wxpm = amrex::Real(-2.0) * (vel(i-1,j-1,k,2)) * idx;
            wxmp = amrex::Real(-2.0) * (vel(i-1,j,k-1,2)) * idx;
            wxmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,2)) * idx;
        }
        else if (bc_type[0][1] == 1) { // slip
            uxpp = amrex::Real(-2.0) * (vel(i-1,j,k,0)) * idx;
            uxpm = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idx;
            uxmp = amrex::Real(-2.0) * (vel(i-1,j,k-1,0)) * idx;
            uxmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,0)) * idx;
            vxpp = amrex::Real(0.0);
            vxpm = amrex::Real(0.0);
            vxmp = amrex::Real(0.0);
            vxmm = amrex::Real(0.0);
            wxpp = amrex::Real(0.0);
            wxpm = amrex::Real(0.0);
            wxmp = amrex::Real(0.0);
            wxmm = amrex::Real(0.0);
        }
    }

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

    // y-wall (lo)
    if (j==domlo.y) {
        if (bc_type[1][0] == 2) { // no-slip
            uypp = amrex::Real(2.0) * (vel(i,j,k,0)) * idy;
            uypm = amrex::Real(2.0) * (vel(i-1,j,k,0)) * idy;
            uymp = amrex::Real(2.0) * (vel(i,j,k-1,0)) * idy;
            uymm = amrex::Real(2.0) * (vel(i-1,j,k-1,0)) * idy;
            vypp = amrex::Real(2.0) * (vel(i,j,k,1)) * idy;
            vypm = amrex::Real(2.0) * (vel(i-1,j,k,1)) * idy;
            vymp = amrex::Real(2.0) * (vel(i,j,k-1,1)) * idy;
            vymm = amrex::Real(2.0) * (vel(i-1,j,k-1,1)) * idy;
            wypp = amrex::Real(2.0) * (vel(i,j,k,2)) * idy;
            wypm = amrex::Real(2.0) * (vel(i-1,j,k,2)) * idy;
            wymp = amrex::Real(2.0) * (vel(i,j,k-1,2)) * idy;
            wymm = amrex::Real(2.0) * (vel(i-1,j,k-1,2)) * idy;
        }
        else if (bc_type[1][0] == 1) { // slip
            uypp = amrex::Real(0.0);
            uypm = amrex::Real(0.0);
            uymp = amrex::Real(0.0);
            uymm = amrex::Real(0.0);
            vypp = amrex::Real(2.0) * (vel(i,j,k,1)) * idy;
            vypm = amrex::Real(2.0) * (vel(i-1,j,k,1)) * idy;
            vymp = amrex::Real(2.0) * (vel(i,j,k-1,1)) * idy;
            vymm = amrex::Real(2.0) * (vel(i-1,j,k-1,1)) * idy;
            wypp = amrex::Real(0.0);
            wypm = amrex::Real(0.0);
            wymp = amrex::Real(0.0);
            wymm = amrex::Real(0.0);
        }
    }
    // y-wall (hi)
    if (j==domhi.y + 1) {
        if (bc_type[1][1] == 2) { // no-slip
            uypp = amrex::Real(-2.0) * (vel(i,j-1,k,0)) * idy;
            uypm = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idy;
            uymp = amrex::Real(-2.0) * (vel(i,j-1,k-1,0)) * idy;
            uymm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,0)) * idy;
            vypp = amrex::Real(-2.0) * (vel(i,j-1,k,1)) * idy;
            vypm = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idy;
            vymp = amrex::Real(-2.0) * (vel(i,j-1,k-1,1)) * idy;
            vymm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,1)) * idy;
            wypp = amrex::Real(-2.0) * (vel(i,j-1,k,2)) * idy;
            wypm = amrex::Real(-2.0) * (vel(i-1,j-1,k,2)) * idy;
            wymp = amrex::Real(-2.0) * (vel(i,j-1,k-1,2)) * idy;
            wymm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,2)) * idy;
        }
        else if (bc_type[1][1] == 1) { // slip
            uypp = amrex::Real(0.0);
            uypm = amrex::Real(0.0);
            uymp = amrex::Real(0.0);
            uymm = amrex::Real(0.0);
            vypp = amrex::Real(-2.0) * (vel(i,j-1,k,1)) * idy;
            vypm = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idy;
            vymp = amrex::Real(-2.0) * (vel(i,j-1,k-1,1)) * idy;
            vymm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,1)) * idy;
            wypp = amrex::Real(0.0);
            wypm = amrex::Real(0.0);
            wymp = amrex::Real(0.0);
            wymm = amrex::Real(0.0);
        }
    }

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

    // z-wall (lo)
    if (k==domlo.z) {
        if (bc_type[2][0] == 2) { // no-slip
            uzpp = amrex::Real(2.0) * (vel(i,j,k,0)) * idz;
            uzpm = amrex::Real(2.0) * (vel(i-1,j,k,0)) * idz;
            uzmp = amrex::Real(2.0) * (vel(i,j-1,k,0)) * idz;
            uzmm = amrex::Real(2.0) * (vel(i-1,j-1,k,0)) * idz;
            vzpp = amrex::Real(2.0) * (vel(i,j,k,1)) * idz;
            vzpm = amrex::Real(2.0) * (vel(i-1,j,k,1)) * idz;
            vzmp = amrex::Real(2.0) * (vel(i,j-1,k,1)) * idz;
            vzmm = amrex::Real(2.0) * (vel(i-1,j-1,k,1)) * idz;
            wzpp = amrex::Real(2.0) * (vel(i,j,k,2)) * idz;
            wzpm = amrex::Real(2.0) * (vel(i-1,j,k,2)) * idz;
            wzmp = amrex::Real(2.0) * (vel(i,j-1,k,2)) * idz;
            wzmm = amrex::Real(2.0) * (vel(i-1,j-1,k,2)) * idz;
        }
        else if (bc_type[2][0] == 1) { // slip
            uzpp = amrex::Real(0.0);
            uzpm = amrex::Real(0.0);
            uzmp = amrex::Real(0.0);
            uzmm = amrex::Real(0.0);
            vzpp = amrex::Real(0.0);
            vzpm = amrex::Real(0.0);
            vzmp = amrex::Real(0.0);
            vzmm = amrex::Real(0.0);
            wzpp = amrex::Real(2.0) * (vel(i,j,k,2)) * idz;
            wzpm = amrex::Real(2.0) * (vel(i-1,j,k,2)) * idz;
            wzmp = amrex::Real(2.0) * (vel(i,j-1,k,2)) * idz;
            wzmm = amrex::Real(2.0) * (vel(i-1,j-1,k,2)) * idz;
        }
    }
    // z-wall (hi)
    if (k==domhi.z + 1) {
        if (bc_type[2][1] == 2) { // no-slip
            uzpp = amrex::Real(-2.0) * (vel(i,j,k-1,0)) * idz;
            uzpm = amrex::Real(-2.0) * (vel(i-1,j,k-1,0)) * idz;
            uzmp = amrex::Real(-2.0) * (vel(i,j-1,k-1,0)) * idz;
            uzmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,0)) * idz;
            vzpp = amrex::Real(-2.0) * (vel(i,j,k-1,1)) * idz;
            vzpm = amrex::Real(-2.0) * (vel(i-1,j,k-1,1)) * idz;
            vzmp = amrex::Real(-2.0) * (vel(i,j-1,k-1,1)) * idz;
            vzmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,1)) * idz;
            wzpp = amrex::Real(-2.0) * (vel(i,j,k-1,2)) * idz;
            wzpm = amrex::Real(-2.0) * (vel(i-1,j,k-1,2)) * idz;
            wzmp = amrex::Real(-2.0) * (vel(i,j-1,k-1,2)) * idz;
            wzmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,2)) * idz;
        }
        else if (bc_type[2][1] == 1) { // slip
            uzpp = amrex::Real(0.0);
            uzpm = amrex::Real(0.0);
            uzmp = amrex::Real(0.0);
            uzmm = amrex::Real(0.0);
            vzpp = amrex::Real(0.0);
            vzpm = amrex::Real(0.0);
            vzmp = amrex::Real(0.0);
            vzmm = amrex::Real(0.0);
            wzpp = amrex::Real(-2.0) * (vel(i,j,k-1,2)) * idz;
            wzpm = amrex::Real(-2.0) * (vel(i-1,j,k-1,2)) * idz;
            wzmp = amrex::Real(-2.0) * (vel(i,j-1,k-1,2)) * idz;
            wzmm = amrex::Real(-2.0) * (vel(i-1,j-1,k-1,2)) * idz;
        }
    }

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

#elif (AMREX_SPACEDIM == 2)
    
    Real uxp = (vel(i,j,k,0) - vel(i-1,j,k,0)) * idx;
    Real uxm = (vel(i,j-1,k,0) - vel(i-1,j-1,k,0)) * idx;
    Real vxp = (vel(i,j,k,1) - vel(i-1,j,k,1)) * idx;
    Real vxm = (vel(i,j-1,k,1) - vel(i-1,j-1,k,1)) * idx;

    // x-wall (lo)
    if (i==domlo.x) {
        if (bc_type[0][0] == 2) { // no-slip
            uxp = amrex::Real(2.0) * (vel(i,j,k,0)) * idx;
            uxm = amrex::Real(2.0) * (vel(i,j-1,k,0)) * idx;
            vxp = amrex::Real(2.0) * (vel(i,j,k,1)) * idx;
            vxm = amrex::Real(2.0) * (vel(i,j-1,k,1)) * idx;
        }
        else if (bc_type[0][0] == 1) { // slip
            uxp = amrex::Real(2.0) * (vel(i,j,k,0)) * idx;
            uxm = amrex::Real(2.0) * (vel(i,j-1,k,0)) * idx;
            vxp = amrex::Real(0.0);
            vxm = amrex::Real(0.0);
        }
    }
    // x-wall (hi)
    if (i==domhi.x + 1) {
        if (bc_type[0][1] == 2) { // no-slip
            uxp = amrex::Real(-2.0) * (vel(i-1,j,k,0)) * idx;
            uxm = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idx;
            vxp = amrex::Real(-2.0) * (vel(i-1,j,k,1)) * idx;
            vxm = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idx;
        }
        else if (bc_type[0][1] == 1) { // slip
            uxp = amrex::Real(-2.0) * (vel(i-1,j,k,0)) * idx;
            uxm = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idx;
            vxp = amrex::Real(0.0);
            vxm = amrex::Real(0.0);
        }
    }

    Real uyp = (vel(i,j,k,0) - vel(i,j-1,k,0)) * idy;
    Real uym = (vel(i-1,j,k,0) - vel(i-1,j-1,k,0)) * idy;
    Real vyp = (vel(i,j,k,1) - vel(i,j-1,k,1)) * idy;
    Real vym = (vel(i-1,j,k,1) - vel(i-1,j-1,k,1)) * idy;

    // y-wall (lo)
    if (j==domlo.y) {
        if (bc_type[1][0] == 2) { // no-slip
            uyp = amrex::Real(2.0) * (vel(i,j,k,0)) * idy;
            uym = amrex::Real(2.0) * (vel(i-1,j,k,0)) * idy;
            vyp = amrex::Real(2.0) * (vel(i,j,k,1)) * idy;
            vym = amrex::Real(2.0) * (vel(i-1,j,k,1)) * idy;
        }
        else if (bc_type[1][0] == 1) { // slip
            uyp = amrex::Real(0.0);
            uym = amrex::Real(0.0);
            vyp = amrex::Real(2.0) * (vel(i,j,k,1)) * idy;
            vym = amrex::Real(2.0) * (vel(i-1,j,k,1)) * idy;
        }
    }
    // y-wall (hi)
    if (j==domhi.y + 1) {
        if (bc_type[1][1] == 2) { // no-slip
            uyp = amrex::Real(-2.0) * (vel(i,j-1,k,0)) * idy;
            uym = amrex::Real(-2.0) * (vel(i-1,j-1,k,0)) * idy;
            vyp = amrex::Real(-2.0) * (vel(i,j-1,k,1)) * idy;
            vym = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idy;
        }
        else if (bc_type[1][1] == 1) { // slip
            uyp = amrex::Real(0.0);
            uym = amrex::Real(0.0);
            vyp = amrex::Real(-2.0) * (vel(i,j-1,k,1)) * idy;
            vym = amrex::Real(-2.0) * (vel(i-1,j-1,k,1)) * idy;
        }
    }
    
    Real ux = amrex::Real(0.5)*(uxp + uxm);
    Real uy = amrex::Real(0.5)*(uyp + uym);
    Real vx = amrex::Real(0.5)*(vxp + vxm);
    Real vy = amrex::Real(0.5)*(vyp + vym);

    return std::sqrt(amrex::Real(2.0) * ux*ux + amrex::Real(2.0) * vy*vy
                     + (uy+vx)*(uy+vx));
#endif

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_cc_to_nodal (int i, int j, int k, amrex::Array4<amrex::Real const> const& rho, bool use_harmonic) noexcept
{
    using namespace amrex;

    Real nodal = 0.0;

#if (AMREX_SPACEDIM == 3)

    if (!use_harmonic) {

        nodal = Real(0.125)*(  rho(i-1,j-1,k-1) + rho(i,j-1,k-1)  + 
                               rho(i-1,j,k-1)   + rho(i,j,k-1)    +
                               rho(i-1,j-1,k)   + rho(i,j-1,k)    +
                               rho(i-1,j,k)     + rho(i,j,k)      );

    }

    else {

        nodal = Real(8.0) / (  (amrex::Real(1.0)/rho(i-1,j-1,k-1)) + (amrex::Real(1.0)/rho(i,j-1,k-1))  + 
                               (amrex::Real(1.0)/rho(i-1,j,k-1))   + (amrex::Real(1.0)/rho(i,j,k-1))    +
                               (amrex::Real(1.0)/rho(i-1,j-1,k))   + (amrex::Real(1.0)/rho(i,j-1,k))    +
                               (amrex::Real(1.0)/rho(i-1,j,k))     + (amrex::Real(1.0)/rho(i,j,k))      );
    }
    return nodal;
#elif (AMREX_SPACEDIM == 2)
    
    if (!use_harmonic) {

        nodal = Real(0.25)*(   rho(i-1,j-1,k)   + rho(i,j-1,k)    +
                               rho(i-1,j,k)     + rho(i,j,k)      );
    }

    else {

        nodal = Real(4.0)/(   (amrex::Real(1.0)/rho(i-1,j-1,k))   + (amrex::Real(1.0)/rho(i,j-1,k))    +
                              (amrex::Real(1.0)/rho(i-1,j,k))     + (amrex::Real(1.0)/rho(i,j,k))      );

    }
    return nodal;
#endif
}

// get nodal concentration from cell centered densities
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_nodal_conc (int i, int j, int k, amrex::Array4<amrex::Real const> const& rho, amrex::Real rho0, amrex::Real rho1) noexcept
{
    using namespace amrex;
    
    Real nodal_conc;

#if (AMREX_SPACEDIM == 3)
    
    Real c_ppp = ((rho1/rho(i,j,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_ppm = ((rho1/rho(i,j,k-1)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_pmp = ((rho1/rho(i,j-1,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_pmm = ((rho1/rho(i,j-1,k-1)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mpp = ((rho1/rho(i-1,j,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mpm = ((rho1/rho(i-1,j,k-1)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mmp = ((rho1/rho(i-1,j-1,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mmm = ((rho1/rho(i-1,j-1,k-1)) - 1.0)/((rho1/rho0) - 1.0);

    nodal_conc = Real(0.125)*(c_ppp+c_ppm+c_pmp+c_pmm+c_mpp+c_mpm+c_mmp+c_mmm);

    return nodal_conc;

#elif (AMREX_SPACEDIM == 2)

    Real c_pp = ((rho1/rho(i,j,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_pm = ((rho1/rho(i,j-1,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mp = ((rho1/rho(i-1,j,k)) - 1.0)/((rho1/rho0) - 1.0);
    Real c_mm = ((rho1/rho(i-1,j-1,k)) - 1.0)/((rho1/rho0) - 1.0);

    nodal_conc = Real(0.25)*(c_pp+c_pm+c_mp+c_mm);

    return nodal_conc;
    
#endif
}

// average cc to edge parallel to y (in 2D) and z in 3D)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_cc_to_edge (int i, int j, int k, amrex::Array4<amrex::Real const> const& rho) noexcept
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 3)

    return Real(0.25)*(  rho(i-1,j-1,k)   + rho(i,j-1,k)    +
                         rho(i-1,j,k)     + rho(i,j,k) );
#elif (AMREX_SPACEDIM == 2)
    
    return Real(0.5)*(   rho(i-1,j,k)   + rho(i,j,k)  );

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

//AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//amrex::Real incflo_get_hydrostatic_pressure(int i, int j, int k, amrex::Array4<amrex::Real const> const& rho,
//                                     amrex::Real /*top*/, int n_top, amrex::Real dx, amrex::Real grav, amrex::Real amb_press) noexcept
//{
//    using namespace amrex;
//    Real hyd_press = amb_press;
//#if (AMREX_SPACEDIM == 3)
//    if (k == n_top) {
//        return hyd_press;
//    }
//    else {
//        for (int kk = n_top-1; kk >= k; --kk) {
//            Real dens = incflo_cc_to_edge(i,j,kk,rho);
//            hyd_press -= dens*grav*dx;
//        }
//        return hyd_press;
//    }
//#elif (AMREX_SPACEDIM == 2)
//    if (j == n_top) {
//        return hyd_press;
//    }
//    else {
//        for (int jj = n_top-1; jj >= j; --jj) {
//            Real dens = incflo_cc_to_edge(i,jj,k,rho);
//            hyd_press -= dens*grav*dx;
//        }
//        return hyd_press;
//    }
//#endif
//}
