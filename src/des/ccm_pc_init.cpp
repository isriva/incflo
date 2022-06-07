#include "ccm_pc.H"

using namespace amrex;

void CCMParticleContainer::InitParticlesAscii (const std::string& particle_init_type, const std::string& file)
{
  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())
  {
    std::ifstream ifs;
    ifs.open(file.c_str(), std::ios::in);

    if (!ifs.good())
      amrex::FileOpenFailed(file);

    int nc = -1;
    ifs >> nc >> std::ws;

    Print() << "Number of cells to be read by input file: " << nc << "\n";

    // Issue an error if nparticles = 0 is specified
    if ( nc == -1 ){
      Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                   \nPerhaps you forgot to specify the number of particles on the first line??? ");
    }

    // we add all the particles to grid 0 and tile 0 and let
    // Redistribute() put them in the right places.
    const int lev  = 0;
    const int grid = 0;
    const int tile = 0;

    auto& particles = DefineAndReturnParticleTile(lev,grid,tile);
    particles.resize(nc);

    Gpu::HostVector<ParticleType> particles_in(nc);

    for (int i = 0; i < nc; i++)
    {
      // We will define the location of the cell to be its centroid
      Real cent_x(0.);
      Real cent_y(0.);

      if (particle_init_type == "particle_locations")
      {
          // amrex::Print() << "... Reading " << NPC << " particles from the inputs file " << std::endl;
          for (int j = 0; j < NPC; j++)
          {
              Real x_0, y_0;
              ifs >> x_0;
              ifs >> y_0;

              if (reverse_x)
                  x_0 *= -1;
              if (reverse_y)
                  y_0 *= -1;

              particles_in[i].rdata(realData::pos_x0+     j *AMREX_SPACEDIM) = x_0;
              particles_in[i].rdata(realData::pos_y0+     j *AMREX_SPACEDIM) = y_0;

              cent_x += x_0;
              cent_y += y_0;
          }

          cent_x /= NPC;
          cent_y /= NPC;

          particles_in[i].pos(0) = cent_x;
          particles_in[i].pos(1) = cent_y;

          for (int j = 0; j < NPC; j++)
          {
              particles_in[i].rdata(realData::pos_x0+     j *AMREX_SPACEDIM) -= cent_x;
              particles_in[i].rdata(realData::pos_y0+     j *AMREX_SPACEDIM) -= cent_y;
          }

      } else {
         amrex::Abort("Dont know this particle_init_type");
      }

      // Here we hard-wire the initial velocity to zero (for convenience)
      for (int j = 0; j < NPC; j++)
      {

          particles_in[i].rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM) = 0.0;
          particles_in[i].rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM) = 0.0;
      }

      // Set id and cpu for this particle
      particles_in[i].id()  = ParticleType::NextID();
      particles_in[i].cpu() = ParallelDescriptor::MyProc();

      // Set other particle properties
      particles_in[i].idata(intData::ccell)     = i;

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }

    auto& aos = particles.GetArrayOfStructs();
    Gpu::DeviceVector<ParticleType>& gpu_particles = aos();

    // Copy particles from host to device
    Gpu::copyAsync(Gpu::hostToDevice, particles_in.begin(), particles_in.end(), gpu_particles.begin());
  }

  Redistribute();
}

