#ifndef INCFLO_STRUCT_FACT_H_
#define INCFLO_STRUCT_FACT_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class IncfloStructFact
{
public:
    IncfloStructFact () = default;

    void define(amrex::BoxArray const& ba,
                amrex::DistributionMapping const& dm,
                amrex::Geometry const& geom,
                amrex::Vector<std::string> const& var_names,
                amrex::Vector<amrex::Real> const& var_scaling,
                int verbosity = 0);

    void define(amrex::BoxArray const& ba,
                amrex::DistributionMapping const& dm,
                amrex::Geometry const& geom,
                amrex::Vector<std::string> const& var_names,
                amrex::Vector<amrex::Real> const& var_scaling,
                amrex::Vector<int> const& pair_a,
                amrex::Vector<int> const& pair_b,
                int verbosity = 0);

    void sample(amrex::MultiFab const& variables, int reset = 0);
    void reset();

    void writePlotFile(int step, amrex::Real time,
                       std::string const& plotfile_base,
                       int zero_mode = 1);

    void finalize(amrex::MultiFab& cov_real_in,
                  amrex::MultiFab& cov_imag_in,
                  int zero_mode = 1);
    void callFinalize(int zero_mode = 1);
    void integrateShells(int step, std::string const& output_base) const;
    void integrateTensorShells(int step, std::string const& output_base) const;

    void writeCheckpoint(std::string const& checkpoint_dir) const;
    void readCheckpoint(std::string const& checkpoint_dir);

    [[nodiscard]] int samples() const noexcept { return m_nsamples; }
    [[nodiscard]] int nCovariances() const noexcept { return m_ncov; }
    [[nodiscard]] amrex::Vector<std::string> const& names() const noexcept { return m_cov_names; }

    void computeFFT(amrex::MultiFab const& variables,
                    amrex::MultiFab& variables_dft_real,
                    amrex::MultiFab& variables_dft_imag,
                    bool unpack = true) const;
    void shiftFFT(amrex::MultiFab& dft_out, int zero_mode = 1) const;

private:
    int m_nvar = 1;
    int m_nvar_unique = 1;
    int m_ncov = 1;
    int m_verbosity = 0;
    int m_nsamples = 0;

    amrex::Vector<amrex::Real> m_scaling;
    amrex::Vector<std::string> m_cov_names;
    amrex::Vector<int> m_pair_a;
    amrex::Vector<int> m_pair_b;
    amrex::Vector<int> m_var_unique;

    amrex::Geometry m_kspace_geom;
    amrex::IntVect m_ncell{AMREX_D_DECL(1, 1, 1)};

    amrex::MultiFab m_cov_real;
    amrex::MultiFab m_cov_imag;
    amrex::MultiFab m_cov_mag;
};

#endif
