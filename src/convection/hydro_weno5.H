#ifndef HYDRO_WENO5_H_
#define HYDRO_WENO5_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

namespace WENO5 {

void ComputeEdgeState(amrex::Box const& bx,
                      AMREX_D_DECL(amrex::Array4<amrex::Real> const& xedge,
                                   amrex::Array4<amrex::Real> const& yedge,
                                   amrex::Array4<amrex::Real> const& zedge),
                      amrex::Array4<amrex::Real const> const& q,
                      int ncomp,
                      AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                   amrex::Array4<amrex::Real const> const& vmac,
                                   amrex::Array4<amrex::Real const> const& wmac),
                      amrex::Box const& domain,
                      amrex::Vector<amrex::BCRec> const& bcs,
                      amrex::BCRec const* d_bcrec,
                      bool is_velocity,
                      amrex::Array4<int const> const& bc_arr = {});

void ExtrapVelToFaces(amrex::MultiFab const& vel,
                      AMREX_D_DECL(amrex::MultiFab& umac,
                                   amrex::MultiFab& vmac,
                                   amrex::MultiFab& wmac),
                      amrex::Geometry const& geom,
                      amrex::Vector<amrex::BCRec> const& h_bcrec,
                      amrex::BCRec const* d_bcrec,
                      bool allow_inflow_on_outflow = false);

void ExtrapVelToFacesBox(AMREX_D_DECL(amrex::Box const& ubx,
                                      amrex::Box const& vbx,
                                      amrex::Box const& wbx),
                         AMREX_D_DECL(amrex::Array4<amrex::Real> const& u,
                                      amrex::Array4<amrex::Real> const& v,
                                      amrex::Array4<amrex::Real> const& w),
                         amrex::Array4<amrex::Real const> const& vcc,
                         amrex::Geometry const& geom,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec,
                         bool allow_inflow_on_outflow = false);

}

#endif
