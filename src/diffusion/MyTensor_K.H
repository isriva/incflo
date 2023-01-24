#ifndef MYTENSOR_K_H_
#define MYTENSOR_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BndryData.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi) noexcept
{
    return (vel(i,j+1,k,n)+vel(i-1,j+1,k,n)-vel(i,j-1,k,n)-vel(i-1,j-1,k,n))*(Real(0.25)*dyi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi) noexcept
{
    return (vel(i+1,j,k,n)+vel(i+1,j-1,k,n)-vel(i-1,j,k,n)-vel(i-1,j-1,k,n))*(Real(0.25)*dxi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi) noexcept
{
    return (vel(i,j,k,n)-vel(i-1,j,k,n))*(Real(1.0)*dxi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi) noexcept
{
    return (vel(i,j,k,n)-vel(i,j-1,k,n))*(Real(1.0)*dyi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi,
                           Array4<Real const> const& bvxlo, Array4<Real const> const& bvxhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddx;
    if (i == dlo.x) {
        if (bct(Orientation::xlo(),n) == AMREX_LO_DIRICHLET && bvxlo) {
            ddx = (bvxlo(i,j,k,n)-bvxlo(i-1,j,k,n))*(Real(0.5)*dxi); //IS need to fix - extrapolate
        } else if (bct(Orientation::xlo(),n) == AMREX_LO_NEUMANN) {
            ddx = Real(0.);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(2.0)*vel(i,j,k,n)*dxi;
        }
    } else if (i == dhi.x+1) {
        if (bct(Orientation::xhi(),n) == AMREX_LO_DIRICHLET && bvxhi) {
            ddx = (bvxhi(i+1,j,k,n)-bvxhi(i,j,k,n))*(Real(0.5)*dxi);  //IS need to fix - extrapolate
        } else if (bct(Orientation::xhi(),n) == AMREX_LO_NEUMANN) {
            ddx = Real(0.);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(-2.0)*vel(i-1,j,k,n)*dxi;
        }
    } else {
        ddx = mltensor_dx_on_xface(i,j,k,n,vel,dxi);
    }
    return ddx;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi,
                           Array4<Real const> const& bvxlo, Array4<Real const> const& bvxhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddy;
    if (i == dlo.x) {
        if (bct(Orientation::xlo(),n) == AMREX_LO_DIRICHLET && bvxlo) {
            if (j == dlo.y) {
                ddy = (bvxlo(i-1,j  ,k,n) * Real(-1.5) +
                       bvxlo(i-1,j+1,k,n) * Real(2.) +
                       bvxlo(i-1,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvxlo(i-1,j  ,k,n) * Real(-1.5) +
                        bvxlo(i-1,j-1,k,n) * Real(2.) +
                        bvxlo(i-1,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = (bvxlo(i-1,j+1,k,n)-bvxlo(i-1,j-1,k,n))*(Real(0.5)*dyi);
            }
        } else if (bct(Orientation::xlo(),n) == AMREX_LO_NEUMANN) {
            ddy = (vel(i,j+1,k,n)-vel(i,j-1,k,n))*(Real(0.5)*dyi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else if (i == dhi.x+1) {
        if (bct(Orientation::xhi(),n) == AMREX_LO_DIRICHLET && bvxhi) {
            if (j == dlo.y) {
                ddy = (bvxhi(i,j  ,k,n) * Real(-1.5) +
                       bvxhi(i,j+1,k,n) * Real(2.) +
                       bvxhi(i,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvxhi(i,j  ,k,n) * Real(-1.5) +
                        bvxhi(i,j-1,k,n) * Real(2.) +
                        bvxhi(i,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = (bvxhi(i,j+1,k,n)-bvxhi(i,j-1,k,n))*(Real(0.5)*dyi);
            }
        } else if (bct(Orientation::xhi(),n) == AMREX_LO_NEUMANN) {
            ddy = (vel(i-1,j+1,k,n)-vel(i-1,j-1,k,n))*(Real(0.5)*dyi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else {
        ddy = mltensor_dy_on_xface(i,j,k,n,vel,dyi);
    }
    return ddy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi,
                           Array4<Real const> const& bvylo, Array4<Real const> const& bvyhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddy;
    if (j == dlo.y) {
        if (bct(Orientation::ylo(),n) == AMREX_LO_DIRICHLET && bvylo) {
            ddy = (bvylo(i,j,k,n)-bvylo(i,j-1,k,n))*(Real(0.5)*dyi); //IS need to fix - extrapolate
        } else if (bct(Orientation::ylo(),n) == AMREX_LO_NEUMANN) {
            ddy = Real(0.);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(2.0)*vel(i,j,k,n)*dyi;
        }
    } else if (j == dhi.y+1) {
        if (bct(Orientation::yhi(),n) == AMREX_LO_DIRICHLET && bvyhi) {
            ddy = (bvyhi(i,j+1,k,n)-bvyhi(i,j,k,n))*(Real(0.5)*dyi);  //IS need to fix - extrapolate
        } else if (bct(Orientation::yhi(),n) == AMREX_LO_NEUMANN) {
            ddy = Real(0.);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(-2.0)*vel(i,j-1,k,n)*dyi;
        }
    } else {
        ddy = mltensor_dy_on_yface(i,j,k,n,vel,dyi);
    }
    return ddy;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi,
                           Array4<Real const> const& bvylo, Array4<Real const> const& bvyhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddx;
    if (j == dlo.y) {
        if (bct(Orientation::ylo(),n) == AMREX_LO_DIRICHLET && bvylo) {
            if (i == dlo.x) {
                ddx = (bvylo(i  ,j-1,k,n) * Real(-1.5) +
                       bvylo(i+1,j-1,k,n) * Real(2.) +
                       bvylo(i+2,j-1,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvylo(i  ,j-1,k,n) * Real(-1.5) +
                        bvylo(i-1,j-1,k,n) * Real(2.) +
                        bvylo(i-2,j-1,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = (bvylo(i+1,j-1,k,n)-bvylo(i-1,j-1,k,n))*(Real(0.5)*dxi);
            }
        } else if (bct(Orientation::ylo(),n) == AMREX_LO_NEUMANN) {
            ddx = (vel(i+1,j,k,n)-vel(i-1,j,k,n))*(Real(0.5)*dxi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else if (j == dhi.y+1) {
        if (bct(Orientation::yhi(),n) == AMREX_LO_DIRICHLET && bvyhi) {
            if (i == dlo.x) {
                ddx = (bvyhi(i  ,j,k,n) * Real(-1.5) +
                       bvyhi(i+1,j,k,n) * Real(2.) +
                       bvyhi(i+2,j,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvyhi(i  ,j,k,n) * Real(-1.5) +
                        bvyhi(i-1,j,k,n) * Real(2.) +
                        bvyhi(i-2,j,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = (bvyhi(i+1,j,k,n)-bvyhi(i-1,j,k,n))*(Real(0.5)*dxi);
            }
        } else if (bct(Orientation::yhi(),n) == AMREX_LO_NEUMANN) {
            ddx = (vel(i+1,j-1,k,n)-vel(i-1,j-1,k,n))*(Real(0.5)*dxi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else {
        ddx = mltensor_dx_on_yface(i,j,k,n,vel,dxi);
    }
    return ddx;
}
}

#if (AMREX_SPACEDIM == 1)
#include <MyTensor_3D_K.H>
//#elif (AMREX_SPACEDIM == 2)
//#include <MyTensor_2D_K.H>
//#else
//#include <MyTensor_3D_K.H>
#endif

#endif
