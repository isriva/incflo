#include <ccm_pc.H>
// #include <incflo.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

Real CCMParticleContainer::ke;
Real CCMParticleContainer::gamma;
Real CCMParticleContainer::fnb;
Real CCMParticleContainer::kb;
Real CCMParticleContainer::ka;
Real CCMParticleContainer::k_adhesive;
Real CCMParticleContainer::orig_domain_area;
Real CCMParticleContainer::cell_a0;
Real CCMParticleContainer::cell_r0;
Real CCMParticleContainer::constant_r0;
Real CCMParticleContainer::cell_x0;
Real CCMParticleContainer::cell_rmax;
Real CCMParticleContainer::lj_length_scale;
Real CCMParticleContainer::lj_repulsive_cutoff;
Real CCMParticleContainer::lj_attractive_cutoff;

int  CCMParticleContainer::move_cell;
Real CCMParticleContainer::x_const_force;
Real CCMParticleContainer::y_const_force;

bool CCMParticleContainer::reverse_x;
bool CCMParticleContainer::reverse_y;

CCMParticleContainer::CCMParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count>
      (amr_core->GetParGDB(), 2)
{
    ReadStaticParameters();

    this->SetVerbose(0);

    nlev = amr_core->maxLevel() + 1;
}

void CCMParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void CCMParticleContainer::PrintParticleCounts ()
{
  const int lev = 0;
  amrex::AllPrintToFile("load_balance") << "Particles on each box: \n";
  long local_count = 0;
  for (CCMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      long np = pti.numParticles();
      local_count += np;
      amrex::AllPrintToFile("load_balance") << "Box:" << pti.index() << ", count: " << np << std::endl;
    }
  amrex::AllPrintToFile("load_balance") << "Total for this process: " << local_count << std::endl << std::endl;
}


void CCMParticleContainer::printParticles ()
{
    const int lev = 0;
    const auto& plevel = GetParticles(lev);

    for (const auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();

       for (int i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "ccell        = " << particles[i].idata(intData::ccell) << " " << std::endl;
          std::cout << "X of cell centroid = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y of cell centroid = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Locations of vertices:" << std::endl;

          for (int j = 0; j < NPC; j++)
          {
            std::cout << "particle x / y" << j << "  = " << 
                         particles[i].rdata(realData::pos_x0+(j-1)*AMREX_SPACEDIM) << " " << 
                         particles[i].rdata(realData::pos_y0+(j-1)*AMREX_SPACEDIM) << std::endl;
          }

          std::cout << std::endl;
       }
    }
}

void CCMParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
   
    ParmParse pp("ccm");
    pp.get("ke",ke);
    pp.get("gamma",gamma);
    pp.get("fnb",fnb);
    pp.get("kb",kb);
    pp.get("ka",ka);

    move_cell  = -1; 
    pp.query("move_cell" , move_cell);

    x_const_force = 0.0;
    y_const_force = 0.0;
    pp.query("x_const_force",x_const_force);
    pp.query("y_const_force",y_const_force);

    reverse_x = false;
    reverse_y = false;
    pp.query("reverse_x",reverse_x);
    pp.query("reverse_y",reverse_y);

    pp.get("k_adhesive",k_adhesive);

    // This input is optional -- we only need it if the current domain area is  
    // grown relative to the original area which determined a0
    orig_domain_area = -1.;
    pp.query("orig_domain_area",orig_domain_area);
}

void CCMParticleContainer::InitParticleParams ()
{
    // We have to read in the particles before we can define a0, x0 and r0
    // cell_a0 = target area of individal cell
    // cell_x0 = distance between vertices of regular polygon with NPC vertices and area cell_a0
    // cell_r0 = distance from centroid    of regular polygon with NPC vertices and area cell_a0
    // constant_r0 = distance from centroid of regular polygon independent of NPC vertices and area cell_a0
    long total_num = TotalNumberOfParticles();

    // If we haven't set orig_domain_area from the inputs file then define it here;
    //    otherwise use the input value
    if (orig_domain_area < 0.)
       orig_domain_area = Geom(0).ProbSize();

    cell_a0 = orig_domain_area / total_num;

    // Formula for area of convex polygon: cell_a0 = 1/2 n R^2 sin (2 pi / n) where n = NPC 
    // so R = sqrt (2 * cell_a0 / n / sin (2 pi / n))
    Real tpi = 2.0 * M_PI;
    Real NN = static_cast<Real>(NPC);
    cell_r0 = std::sqrt(2.0 * cell_a0 / NN / sin (tpi / NN) );

    // R that stays constant regardless of NPC
    constant_r0 = std::sqrt(cell_a0 / M_PI);
    amrex::Print() << "constant_r0:  " << constant_r0 << std::endl;

    // Two equivalent formulae for the side length x_0
    cell_x0 = 2.0 * cell_r0     * std::sin(M_PI/NN);
    cell_x0 = 2.0 * constant_r0 * std::sqrt(M_PI/NN * std::tan(M_PI/NN));

    if ( std::abs(2.0 * cell_r0     * std::sin(M_PI/NN)) - 2.0 * constant_r0 * std::sqrt(M_PI/NN * std::tan(M_PI/NN)) 
         > 1.e-8 )
    {
       amrex::Print() << "Cell_x0 defined as " << 2.0 * cell_r0     * std::sin(M_PI/NN) << " but also as " <<
                          2.0 * constant_r0 * std::sqrt(M_PI/NN * std::tan(M_PI/NN)) << std::endl;
       amrex::Abort("Definitions of cell_x0 dont match!");
    }

    // Maximum distance from a cell centroid to a cell vertex -- over all vertices and all cells
    cell_rmax = cell_r0;

    lj_length_scale      = 0.15578876 * constant_r0;
    lj_repulsive_cutoff  = 0.07789438 * constant_r0;
    lj_attractive_cutoff = 0.62315503 * constant_r0;

    amrex::Print() << " " << std::endl;
    amrex::Print() << "************************************************************ " << std::endl;
    amrex::Print() << "Number of particles read in:                                 " << total_num << std::endl;
    if (orig_domain_area == Geom(0).ProbSize())
         amrex::Print() << "Total domain area:                                           " << orig_domain_area << std::endl;
    else 
    {
         amrex::Print() << "Total area of original domain:                               " << orig_domain_area << std::endl;
         amrex::Print() << "Total area of current  domain:                               " << Geom(0).ProbSize() << std::endl;
    }
    amrex::Print() << "Target area of each cell:                                    " << cell_a0 << std::endl;
    amrex::Print() << "Distance from centroid    of regular polygon with this area: " << cell_r0 << std::endl;
    amrex::Print() << "Distance between vertices of regular polygon with this area: " << cell_x0 << std::endl;
    amrex::Print() << "Cutoff of repulsive  Leonard-Jones force                     " << lj_repulsive_cutoff << std::endl;
    amrex::Print() << "Cutoff of attractive Leonard-Jones force                     " << lj_attractive_cutoff << std::endl;
    amrex::Print() << "************************************************************ " << std::endl;
    amrex::Print() << " " << std::endl;

    // Sanity check
    if (std::abs(cell_a0 - (0.5 * NN * cell_r0 * cell_r0 * sin (tpi / NN))) > 1.e-12 )
    {
       amrex::Print() << "Area formula gives " << (0.5 * NN * cell_r0 * cell_r0 * sin (tpi / NN)) << std::endl;
       amrex::Abort("Ooops -- got area formula wrong");
    }
}
