#include <ccm_pc.H>
// #include <incflo.H>

using namespace amrex;

void CCMParticleContainer::WriteToAscii (const std::string& vtk_filename,int nstep, Real time)
{
  BL_PROFILE("CCMParticleContainer::WriteToAscii");
   
  const std::string& filename = amrex::Concatenate(vtk_filename,nstep,6);

  amrex::Print()  << " "                                              << std::endl;
  amrex::Print()  << "  Writing vtkfile " << filename <<  " at time " << time << std::endl;
  amrex::Print()  << " "                                              << std::endl;

  const int lev = 0;
  const auto& plevel = GetParticles(lev);
  long total_number_of_particles = NumberOfParticlesAtLevel(lev);
  long total_number_of_vertices  = NPC * total_number_of_particles;
    
  // First write all the header information just once (from the I/O procesor)
  if (ParallelDescriptor::IOProcessor())
  {
     //
     // Have I/O processor open file and write out particle metadata.
     //
     std::ofstream File;
     File.open(filename + ".vtk", std::ios::out | std::ios::trunc);
       
     if (!File.good())
        amrex::FileOpenFailed(filename+".vtk");

     // vtkAscii File contents
     File << "# vtk DataFile Version 3.0 "         << std::endl;
     File << "hexagon using polygon function vtk " << std::endl;
     File << "ASCII"                               << std::endl;
     File << "DATASET POLYDATA"                    << std::endl;
     File << " "                                   << std::endl;
     File << "POINTS" << " " << total_number_of_vertices << " " << "float" << std::endl;
     File.flush();
     File.close();
     if (!File.good())
         amrex::Abort("CCMParticleContainer::WriteVTKAscii():  problem writing file");
  }

  ParallelDescriptor::Barrier();

  const int MyProc = ParallelDescriptor::MyProc();

  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".vtk", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".vtk");

          for (const auto& kv : plevel)
          {
            const auto& particles = kv.second.GetArrayOfStructs();

            for (int i = 0; i < particles.numParticles(); ++i)
            {
               for (int j = 0; j < NPC; j++)
               {
                  // To print the absolute position we start with the centroid and add the offset from the centroid
                  Real x = particles[i].pos(0) + particles[i].rdata(realData::pos_x0+(j)*AMREX_SPACEDIM);
                  Real y = particles[i].pos(1) + particles[i].rdata(realData::pos_y0+(j)*AMREX_SPACEDIM);

                  if (reverse_x)
                      x *= -1.;
                  if (reverse_y)
                      y *= -1.;
    
                  File << x << " " << y << " " << 0.00 << std::endl;
               }
               File    << " "                                                            << std::endl;
            }
          } // loop over particles at level "lev"

        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc

    // POLYGONS # of cells #, of vertices per cell   
    if (ParallelDescriptor::IOProcessor())
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream File;
        File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        File.open(filename+".vtk", std::ios::out|std::ios::app);
        if (!File.good())
            amrex::FileOpenFailed(filename+".vtk");

        File << "POLYGONS" << " " <<         total_number_of_particles << " " << 
                                     (NPC+1)*total_number_of_particles << std::endl;

        // Write polygon info
        for (long i = 0; i < total_number_of_particles; ++i)
        {
             File << NPC << " "; 
             for (int j = 0; j < NPC; j++)
                File << i*NPC+j << " "; 
             File     << std::endl;
        }

        File << " "                                                       << std::endl;
        File << "CELL_DATA"     << " "       << total_number_of_particles << std::endl;
        File << "POINT_DATA"    << " "       << total_number_of_vertices  << std::endl; 
        File << "SCALARS "      << "cell "   << "int"                         << std::endl;
        File << "LOOKUP_TABLE " << "default"                                  << std::endl;

        // Write lookup info - this will be used to color the vertices of one cell the same color
        for (long i = 0; i < total_number_of_particles; ++i)
        {
             //File << NPC << " "; 
             for (int j = 0; j < NPC; j++)
                File << i << " "; 
             File     << std::endl;
        }
        File     << std::endl;

        File.flush();
        File.close();
        if (!File.good())
             amrex::Abort("CCMParticleContainer::WriteVTKAscii():  problem writing file");
     } // I/O Proc
     
}

void CCMParticleContainer::WriteForRestart (const std::string& restart_filename, int nstep, Real time)
{
  BL_PROFILE("CCMParticleContainer::WriteForRestart");
   
  const std::string& filename = amrex::Concatenate(restart_filename,nstep,6);

  amrex::Print()  << " "                                                                    << std::endl;
  amrex::Print()  << "  Writing particle restart file " << filename <<  " at time " << time << std::endl;
  amrex::Print()  << " "                                                                    << std::endl;

  const int lev = 0;
  const auto& plevel = GetParticles(lev);
  long total_number_of_particles = NumberOfParticlesAtLevel(lev);
    
  // First write all the header information just once (from the I/O procesor)
  if (ParallelDescriptor::IOProcessor())
  {
     //
     // Have I/O processor open file and write out particle metadata.
     //
     std::ofstream File;
     File.open(filename + ".dat", std::ios::out | std::ios::trunc);
       
     if (!File.good())
        amrex::FileOpenFailed(filename+".dat");

     // File contents
     File << total_number_of_particles << std::endl;
     File.flush();
     File.close();
     if (!File.good())
         amrex::Abort("CCMParticleContainer::WriteVTKAscii():  problem writing file");
  }
  ParallelDescriptor::Barrier();

  const int MyProc = ParallelDescriptor::MyProc();

  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".dat", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".dat");

          for (const auto& kv : plevel)
          {
            const auto& particles = kv.second.GetArrayOfStructs();

            for (int i = 0; i < particles.numParticles(); ++i)
            {
               for (int j = 0; j < NPC; j++)
               {
                  // To print the absolute position we start with the centroid and add the offset from the centroid
                  Real x = particles[i].pos(0) + particles[i].rdata(realData::pos_x0+(j)*AMREX_SPACEDIM);
                  Real y = particles[i].pos(1) + particles[i].rdata(realData::pos_y0+(j)*AMREX_SPACEDIM);

                  // We use these to test the symmetry of our equations
                  if (reverse_x)
                      x *= -1.;
                  if (reverse_y)
                      y *= -1.;
    
                  File << x << " " << y << std::endl;
               }
               File    << " "                                                            << std::endl;
            }
          } // loop over particles at level "lev"

        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc
}
