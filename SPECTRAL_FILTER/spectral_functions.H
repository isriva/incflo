#ifndef SPECTRAL_FUNCTIONS_H_
#define SPECTRAL_FUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>

enum class SpectralFilterType : int {
    Sharp,
    SincSq
};

struct SpectralFilterOptions {
    SpectralFilterType filter_type = SpectralFilterType::Sharp;
    bool zero_outside_range = true;
};

void SpectralReadCheckPoint(const std::string& restart_file,
                            const amrex::Vector<int>& is_periodic,
                            amrex::Geometry& geom,
                            amrex::BoxArray& ba,
                            amrex::DistributionMapping& dmap,
                            amrex::MultiFab& velocity,
                            int& step,
                            amrex::Real& time);

void SpectralVelDecomp(const amrex::MultiFab& velocity,
                       amrex::MultiFab& velocity_filter,
                       amrex::Real kmin,
                       amrex::Real kmax,
                       const SpectralFilterOptions& filter_options,
                       const amrex::Geometry& geom);

void SpectralVelProductDecomp(const amrex::MultiFab& velocity, 
                              amrex::MultiFab& vv_filter, 
                              amrex::Real kmin, 
                              amrex::Real kmax, 
                              const SpectralFilterOptions& filter_options,
                              const amrex::Geometry& geom);

void ProcessDeltaEtaSpectrum(int step, 
                             const amrex::MultiFab& output, 
                             const amrex::Geometry& geom, 
                             amrex::Real kmin, 
                             amrex::Real kmax,
                             const SpectralFilterOptions& filter_options);

void SpectralWritePlotFile(int step,
                           amrex::Real kmin,
                           amrex::Real kmax,
                           const SpectralFilterOptions& filter_options,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& velocity,
                           const amrex::MultiFab& velocity_filter,
                           const amrex::MultiFab& previous_velocity_filter,
                           const amrex::MultiFab& vv_filter,
                           bool use_prime_tau);

void SpectralWriteFourierPlotFile(int step,
                                  amrex::Real kmin,
                                  amrex::Real kmax,
                                  const SpectralFilterOptions& filter_options,
                                  const amrex::Geometry& geom,
                                  const amrex::MultiFab& velocity_filter,
                                  const amrex::MultiFab& vv_filter);

#endif
