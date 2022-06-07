#ifndef _CCM_PC_H_
#define _CCM_PC_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_FabArray.H>
#include <AMReX_FabArray.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#define NPC 20

struct realData
{
  enum { // Real particle data stored in p.m_data
    pos_x0 = 0,
    pos_y0,
    count = 2*AMREX_SPACEDIM*NPC  // Positions and velocity 
  };
};

struct intData
{
  enum { // Integer particle data stored in p.m_data
    ccell = 0,
    count = 1 // Just ccell at this point
  };
};

///
/// We add a particle to the neighbor list if its center is within neighborhood
///
struct CCMCheckPair
{
  amrex::Real neighborhood;
  CCMCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P& p1, const P& p2) const
  { return (p1.pos() - p2.pos()).radSquared() <= neighborhood*neighborhood; }
};

class CCMParticleContainer final

    : public amrex::NeighborParticleContainer<realData::count,intData::count>
{
public:
    using CCMParIter = amrex::ParIter<realData::count,intData::count,0,0>;

    // Get number of particles
    int NumberOfParticles (CCMParIter& pti)
    { return pti.GetArrayOfStructs().size(); }

    CCMParticleContainer (amrex::AmrCore* amr_core);

    virtual ~CCMParticleContainer () {};

    void InitData ();

    void InitParticlesAscii (const std::string& init_type, const std::string& file);

    void InitParticleParams ();

    void WriteToAscii (const std::string& file,int nstep,amrex::Real time);

    void WriteForRestart (const std::string& file,int nstep,amrex::Real time);

    void printParticles ();

    void AllocData ();

    void PrintParticleCounts ();

    void ComputeAndPrintDiagnostics (int lev,
                                    int nstep,
                                    amrex::Real time);

    void EvolveParticles (int lev,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time);

    void CCM_PC_InitCollisionParams ();

    void output (int estatus,
                 int finish,
                 int nstep,
                 amrex::Real dt,
                 amrex::Real time);

    // Parameters inmported from input file
    static int vtk_ascii_int;
    static amrex::Real ke;
    static amrex::Real gamma;
    static amrex::Real fnb;
    static amrex::Real kb;
    static amrex::Real ka;
    static amrex::Real k_adhesive;
    static amrex::Real orig_domain_area;
    static amrex::Real cell_a0;
    static amrex::Real cell_r0;
    static amrex::Real constant_r0;
    static amrex::Real cell_x0;
    static amrex::Real cell_rmax;
    static amrex::Real lj_length_scale;
    static amrex::Real lj_repulsive_cutoff;
    static amrex::Real lj_attractive_cutoff;

    static int move_cell;
    static amrex::Real x_const_force;
    static amrex::Real y_const_force;

    static bool reverse_x;
    static bool reverse_y;

  private:
    int nlev;

    static int domain_bc[6];

    static void ReadStaticParameters ();
};

#endif
