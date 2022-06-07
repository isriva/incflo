#include <ccm_pc.H>
#include <ccm_pc_K.H>

using namespace amrex;

void CCMParticleContainer::EvolveParticles (int lev,
                                            int nstep,
                                            Real subdt,
                                            Real time)
{
    BL_PROFILE_REGION_START("ccm_dem::EvolveParticles()");
    BL_PROFILE("ccm_dem::EvolveParticles()");

    amrex::Print() << "Step " << nstep << ": Evolving particles with dt = " << subdt << std::endl;

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle forces                                            *
     *   -> particle-particle torques                                           *
     ***************************************************************************/
    std::map<PairIndex, Gpu::DeviceVector<Real>> fc;

    std::map<PairIndex, bool> tile_has_walls;

    for (CCMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        fc[index]   = Gpu::DeviceVector<Real>();
    }

    // Within a grid, how far (in physical distance) should we look to consider another cell a neighbor?
    // cell_rmax is the max distance from a centroid to a vertex
    Real max_dist_of_a_neighbor = 2. * cell_rmax + lj_attractive_cutoff;
    
    // Because we defined the CCMParticleContainer with ngrow = 1, we only fill particles in one ghost cell
    // of each patch -- this had better be enough!
    Real dx = this->Geom(0).CellSize()[0];

    amrex::Print() << "dx vs. max interaction distance " 
                   << dx << " " << max_dist_of_a_neighbor << std::endl;

    // Note: the "2" here comes from the fact that we defined
    //CCMParticleContainer::CCMParticleContainer (AmrCore* amr_core)
    //    : NeighborParticleContainer<realData::count,intData::count>
    //      (amr_core->GetParGDB(), 2)

    long total_number_of_particles = NumberOfParticlesAtLevel(0);
    //if (total_number_of_particles  > 1 && max_dist_of_a_neighbor > 2.*dx)
    //  amrex::Abort("dx is too small -- we need more than two ghost cells");

    // How often we re-compute the neighbor lists and redistribute
    int nbor_interval = 10;

    /****************************************************************************
     * Advance by one step in time
     ***************************************************************************/
    {
        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (nstep % nbor_interval == 0) {
            clearNeighbors();
            Redistribute(0, 0, 0, 1);
            fillNeighbors();
            // send in "false" for sort_neighbor_list option

            amrex::Print() << " ... rebuilding neighbor list " << std::endl;
            Real search_radius = max_dist_of_a_neighbor;
            buildNeighborList(CCMCheckPair(search_radius), false);
        } else {
            updateNeighbors();
        }

        // Array2D<int,1,101,1,101> my_matrix;
        //for (int j = 1; j <= 101; j++)
        //   for (int i = 1; i <= 101; i++)
        //       my_matrix(i,j) = 0;

        /********************************************************************
         * Loop over particle tiles at this level                           *
         *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (CCMParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            PairIndex index(pti.index(), pti.LocalTileIndex());

            const int nrp = GetParticles(lev)[index].numRealParticles();

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            // Neighbor particles
#ifdef AMREX_USE_GPU
            int size_ng = aos.numNeighborParticles();
#else
            int size_ng = neighbors[lev][index].size();
#endif

            // Particle-particle forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            fc[index].clear();
            fc[index].resize(nrp*NPC*AMREX_SPACEDIM,0.0);

            Real* fc_ptr = fc[index].dataPtr();

            BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);
            
            // nrp is the number of particles ("particle = cancer cell) 

            /*********************************************************************
             * Particle-Particle spring-dashpot and LJ - SAME CELL ONLY *
             ********************************************************************/
            amrex::ParallelFor(nrp,
                [pstruct,fc_ptr,nrp,subdt]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  ParticleType& p1 = pstruct[i];

                 // Here we iterate over all the vertices of this cell 
                 for (int j = 0; j < NPC; j++)
                 {
                      /********************************************************************
                      * Compute forces on vertex j from vertices (j-1) and (j+1)
                      ********************************************************************/
                      // Location of particle j
                      Real x_0 = p1.rdata(realData::pos_x0+j*AMREX_SPACEDIM);
                      Real y_0 = p1.rdata(realData::pos_y0+j*AMREX_SPACEDIM);
                      Real u_0 = p1.rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM);
                      Real v_0 = p1.rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM);

                      // First compute forces on j from j+1
                      int j_pls = (j == NPC-1) ? 0 : j+1;
                      Real x_pls = p1.rdata(realData::pos_x0+j_pls*AMREX_SPACEDIM);
                      Real y_pls = p1.rdata(realData::pos_y0+j_pls*AMREX_SPACEDIM);
                      Real u_pls = p1.rdata(realData::pos_x0+(NPC+j_pls)*AMREX_SPACEDIM);
                      Real v_pls = p1.rdata(realData::pos_y0+(NPC+j_pls)*AMREX_SPACEDIM);

                      Real fx_pls,fy_pls;
                      computeSpringDashpot(x_0,y_0,u_0,v_0,x_pls,y_pls,u_pls,v_pls,
                                           cell_x0,ke,gamma,constant_r0,fx_pls,fy_pls);

                      // Now compute forces on j from j-1
                      int j_mns = (j == 0) ? NPC-1 : j-1;
                      Real x_mns = p1.rdata(realData::pos_x0+j_mns*AMREX_SPACEDIM);
                      Real y_mns = p1.rdata(realData::pos_y0+j_mns*AMREX_SPACEDIM);
                      Real u_mns = p1.rdata(realData::pos_x0+(NPC+j_mns)*AMREX_SPACEDIM);
                      Real v_mns = p1.rdata(realData::pos_y0+(NPC+j_mns)*AMREX_SPACEDIM);

                      Real fx_mns,fy_mns;
                      computeSpringDashpot(x_0,y_0,u_0,v_0,x_mns,y_mns,u_mns,v_mns,
                                           cell_x0,ke,gamma,constant_r0,fx_mns,fy_mns);

                      // Each particle updates its force (no need for atomics)
                      int force_index = (i*NPC + j)*AMREX_SPACEDIM; 
                      fc_ptr[force_index  ] += fx_mns + fx_pls;
                      fc_ptr[force_index+1] += fy_mns + fy_pls;

                      /********************************************************************
                      * Compute forces on vertex j from line segments that do not include j
                      ********************************************************************/
                      // Loop over all vertices that are *not* j or next to vertex j
                      //    so not equal to j-1,j or j+1
                      for (int j_nb = 0;  j_nb <= NPC - 1; j_nb++)
                      {
                         int j_mns = j_nb;
                         int j_pls = (j_nb == NPC -1) ? 0 : j_nb+1;

                         if (j == j_mns || j == j_pls) continue;

                         Real x_mns = p1.rdata(realData::pos_x0+j_mns*AMREX_SPACEDIM);
                         Real y_mns = p1.rdata(realData::pos_y0+j_mns*AMREX_SPACEDIM);
                         Real x_pls = p1.rdata(realData::pos_x0+j_pls*AMREX_SPACEDIM);
                         Real y_pls = p1.rdata(realData::pos_y0+j_pls*AMREX_SPACEDIM);

                         Real fx,fy,r_dist;
                         computeLeonardJones(x_0,y_0,x_mns,y_mns,x_pls,y_pls,r_dist,
                                             fnb,k_adhesive,fx,fy,
                                             lj_length_scale,lj_repulsive_cutoff,lj_attractive_cutoff,
                                             constant_r0,false);

                         if (fnb > 1.0e10){
                            amrex::Print() << "PARTICLE TO PARTICLE REPULSION FORCE ON PARTICLE " << j << 
                                " FROM " << j_mns << " AND " << j_pls <<  
                                " IN CELL " << i << " IS " << "[" <<fx << "," << fy << "]" << std::endl;
                         }
                         
                         // This is the force on the vertex being pushed away by the line segement
                         int force_index = (i*NPC + j)*AMREX_SPACEDIM; 
                         fc_ptr[force_index  ] += fx;
                         fc_ptr[force_index+1] += fy;

                         // These are the equal-and-opposite forces on the two ends of the line segment
                         int fi_1 = (i*NPC + j_mns)*AMREX_SPACEDIM; 
                         fc_ptr[fi_1  ] -= 0.5*fx;
                         fc_ptr[fi_1+1] -= 0.5*fy;
                         int fi_2 = (i*NPC + j_pls)*AMREX_SPACEDIM; 
                         fc_ptr[fi_2  ] -= 0.5*fx;
                         fc_ptr[fi_2+1] -= 0.5*fy;

                      } // loop over j_nb
                 } // loop over j
              }); // loop over p1

            auto nbor_data = m_neighbor_list[lev][index].data();

            /*********************************************************************
             * Cell to cell repulsion force -- DIFFERENT CELLS ONLY              * 
             * This loop does part of the LJ force that acts on the vertex
             ********************************************************************/

            amrex::ParallelFor(nrp,
                // [pstruct,fc_ptr,nbor_data,nrp,subdt,&my_matrix]
                [pstruct,fc_ptr,nbor_data,nrp,subdt]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  ParticleType& p1 = pstruct[i];
                  int nbor = 0;
                  for (const auto& p2 : nbor_data.getNeighbors(i))
                  {
                      if (p1.idata(intData::ccell) !=  p2.idata(intData::ccell))
                      {
                              
                        //amrex::Print() << "p1.idata(intData::ccell): " << p1.idata(intData::ccell) << std::endl;
                        //amrex::Print() << "p2.idata(intData::ccell): " << p2.idata(intData::ccell) << std::endl;
                        // Here we iterate over all the vertices of this cell 
                        for (int j = 0; j < NPC; j++)
                        {
                           // Location of particle j
                           Real x_0 = p1.rdata(realData::pos_x0+j*AMREX_SPACEDIM) + p1.pos(0);
                           Real y_0 = p1.rdata(realData::pos_y0+j*AMREX_SPACEDIM) + p1.pos(1);

                           for (int j_nb = 0; j_nb <= NPC - 1; j_nb++ )
                           {
                              int j_mns = j_nb;
                              int j_pls = (j_nb == NPC -1) ? 0 : j_nb+1;

                              Real x_mns = p2.rdata(realData::pos_x0+j_mns*AMREX_SPACEDIM) + p2.pos(0);
                              Real y_mns = p2.rdata(realData::pos_y0+j_mns*AMREX_SPACEDIM) + p2.pos(1);
                              Real x_pls = p2.rdata(realData::pos_x0+j_pls*AMREX_SPACEDIM) + p2.pos(0);
                              Real y_pls = p2.rdata(realData::pos_y0+j_pls*AMREX_SPACEDIM) + p2.pos(1);

                              Real fx,fy,r_dist;
                              computeLeonardJones(x_0,y_0,x_mns,y_mns,x_pls,y_pls,r_dist,
                                                  fnb,k_adhesive,fx,fy,
                                                  lj_length_scale,lj_repulsive_cutoff,lj_attractive_cutoff,
                                                  constant_r0,true);

                              //if (std::abs(fx) > 0. || std::abs(fy) > 0.){
                              //if (r_dist < 0.03498328){
                              //   my_matrix(p1.id(),p2.id()) = 1;
                              //}   

                              if (fnb > 1.0e10){
                                 amrex::Print() << "PARTICLE TO PARTICLE REPULSION FORCE ON PARTICLE " << j << 
                                     " IN CELL " << p1.idata(intData::ccell) << " FROM PARTICLES " << j_mns << " AND " << j_pls <<  
                                     " IN CELL " << p2.idata(intData::ccell) << " IS " << "[" <<fx << "," << fy << "]" << std::endl;
                              }
                            
                              int force_index = (i*NPC + j)*AMREX_SPACEDIM; 
                              fc_ptr[force_index  ] += fx;
                              fc_ptr[force_index+1] += fy;
                           } // j_nb in particle 2
                        } // j in particle 1
                        nbor++;
                      } // p1 != p2
                  } // loop over p2
                  // amrex::Print() << "PARTICLE " << i << " HAS " << nbor << " neighbors " << std::endl;
              }); // loop over p1

#if 0
        int cell[102] = {0};
        int bins[12] = {0};
        //int j = 7;
        for (int j = 1; j <= 101; j++){
            for (int i = 1; i <= 101; i++){
                if (my_matrix(i,j) > 0){
                   cell[j] += 1;
                }
                //amrex::Print() << " my_matrix: " << my_matrix(i,j) << std::endl;
            }
            //amrex::Print() << cell[j] << std::endl;
            if (cell[j] == 0)
               bins[0] += 1;
            else if (cell[j] == 1)
               bins[1] += 1;
            else if (cell[j] == 2)
               bins[2] += 1;
            else if (cell[j] == 3)
               bins[3] += 1;
            else if (cell[j] == 4)
               bins[4] += 1;
            else if (cell[j] == 5)
               bins[5] += 1;
            else if (cell[j] == 6)
               bins[6] += 1;
            else if (cell[j] == 7)
               bins[7] += 1;
            else if (cell[j] == 8)
               bins[8] += 1;
            else if (cell[j] == 9)
               bins[9] += 1;
            else if (cell[j] == 10)
               bins[10] += 1;
        }    
        amrex::Print() << " BREAK " << std::endl;
        int binsTotal = 0;
        for (int l = 0; l < 11; l++)
           binsTotal += bins[l];
        amrex::Print() << " 0 sides: " << bins[0] << std::endl;
        amrex::Print() << " 1 sides: " << bins[1] << std::endl;
        amrex::Print() << " 2 sides: " << bins[2] << std::endl;
        amrex::Print() << " 3 sides: " << bins[3] << std::endl;
        amrex::Print() << " 4 sides: " << bins[4] << std::endl;
        amrex::Print() << " 5 sides: " << bins[5] << std::endl;
        amrex::Print() << " 6 sides: " << bins[6] << std::endl;
        amrex::Print() << " 7 sides: " << bins[7] << std::endl;
        amrex::Print() << " 8 sides: " << bins[8] << std::endl;
        amrex::Print() << " 9 sides: " << bins[9] << std::endl;
        amrex::Print() << " 10 sides: " << bins[10] << std::endl;
        amrex::Print() << " binsTotal: " << binsTotal << std::endl;
#endif

            /*********************************************************************
             * Cell to cell repulsion force -- DIFFERENT CELLS ONLY              * 
             * This loop does the equal-and-opposite part of the LJ repulsion force --
             *     the force acting on the line segement 
             ********************************************************************/
            amrex::ParallelFor(nrp,
                [pstruct,fc_ptr,nbor_data,nrp,subdt]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  ParticleType& p1 = pstruct[i];
                  for (const auto& p2 : nbor_data.getNeighbors(i))
                  {
                      if (p1.idata(intData::ccell) !=  p2.idata(intData::ccell))
                      {
                        // Here we iterate over all the line segments of this cell 
                        for (int j = 0; j < NPC; j++)
                        {
                            int j_mns = j;
                            int j_pls = (j == NPC -1) ? 0 : j+1;

                            // Location of end points of line segment between j_mns and j_pls of p1
                            Real x_mns = p1.rdata(realData::pos_x0+j_mns*AMREX_SPACEDIM) + p1.pos(0);
                            Real y_mns = p1.rdata(realData::pos_y0+j_mns*AMREX_SPACEDIM) + p1.pos(1);
                            Real x_pls = p1.rdata(realData::pos_x0+j_pls*AMREX_SPACEDIM) + p1.pos(0);
                            Real y_pls = p1.rdata(realData::pos_y0+j_pls*AMREX_SPACEDIM) + p1.pos(1);

                            // Here we iterate over all the vertices of the neighbor cell
                            for (int j_nb = 0; j_nb <= NPC - 1; j_nb++ )
                            {
                                // Location of vertex j_nb of p2
                                Real x_0 = p2.rdata(realData::pos_x0+j_nb*AMREX_SPACEDIM) + p2.pos(0);
                                Real y_0 = p2.rdata(realData::pos_y0+j_nb*AMREX_SPACEDIM) + p2.pos(1);

                                Real fx,fy,r_dist;
                                computeLeonardJones(x_0,y_0,x_mns,y_mns,x_pls,y_pls,r_dist,
                                                    fnb,k_adhesive,fx,fy,
                                                    lj_length_scale,lj_repulsive_cutoff,lj_attractive_cutoff,
                                                    constant_r0,true);

                                int force_index_m = (i*NPC + j_mns)*AMREX_SPACEDIM; 
                                fc_ptr[force_index_m  ] -= 0.5*fx;
                                fc_ptr[force_index_m+1] -= 0.5*fy;

                                int force_index_p = (i*NPC + j_pls)*AMREX_SPACEDIM; 
                                fc_ptr[force_index_p  ] -= 0.5*fx;
                                fc_ptr[force_index_p+1] -= 0.5*fy;

                            } // j_nb in particle 2
                        } // j_mns,j_pls in particle 1
                      } // p1 != p2
                  } // loop over p2
              }); // loop over p1

            /*********************************************************************
             * Calculate Fb (theta force) -- SAME CELL ONLY
             ********************************************************************/
            amrex::ParallelFor(nrp,
              [pstruct,subdt,fc_ptr]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                ParticleType& p= pstruct[i];
                {
                    Real theta_0[NPC] = {};
                    Real theta_1[NPC] = {};
                    Real x_0_temp[NPC] = {};
                    Real y_0_temp[NPC] = {};
                          
                    // Compute theta_0
                    for (int j = 0; j < NPC; j++)
                    {
                      // Location of particle j
                      Real x_0 = p.rdata(realData::pos_x0+j*AMREX_SPACEDIM);
                      Real y_0 = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);
                      Real u_0 = p.rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM);
                      Real v_0 = p.rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM);
                      
                      // Location of particle j+1
                      int j_pls = (j == NPC-1) ? 0 : j+1;
                      Real x_pls = p.rdata(realData::pos_x0+j_pls*AMREX_SPACEDIM);
                      Real y_pls = p.rdata(realData::pos_y0+j_pls*AMREX_SPACEDIM);

                      int j_mns = (j == 0) ? NPC-1 : j-1;
                      Real x_mns = p.rdata(realData::pos_x0+j_mns*AMREX_SPACEDIM);
                      Real y_mns = p.rdata(realData::pos_y0+j_mns*AMREX_SPACEDIM);

                      Real normal_sum_x, normal_sum_y;

                      computeTheta(x_0, y_0, x_pls, y_pls, x_mns, y_mns, theta_0[j],
                                   normal_sum_x, normal_sum_y);

                      x_0_temp[j] = x_0 + subdt * u_0;
                      y_0_temp[j] = y_0 + subdt * v_0;
                    } 

                    // Compute theta_1
                    for (int j = 0; j < NPC; j++)
                    {
                      
                      Real x_0 = x_0_temp[j];
                      Real y_0 = y_0_temp[j];

                      // Location of particle j+1
                      int j_pls = (j == NPC-1) ? 0 : j+1;
                      Real x_pls = x_0_temp[j_pls];
                      Real y_pls = y_0_temp[j_pls];

                      int j_mns = (j == 0) ? NPC-1 : j-1;
                      Real x_mns = x_0_temp[j_mns];
                      Real y_mns = y_0_temp[j_mns];

                      Real normal_sum_x, normal_sum_y;

                      computeTheta(x_0, y_0, x_pls, y_pls, x_mns, y_mns, theta_1[j],
                                   normal_sum_x, normal_sum_y);

                      Real fb_mag = 0.;
#if 0
                      // Pick a cell to apply a constant force in input file
                      if (p.idata(intData::ccell) == move_cell)
                      {
                         Real kb_small = kb /20;
                         fb_mag = kb_small * tan( (theta_1[j]-theta_0[j]) / 2);

                      } else {

                         fb_mag = kb * tan( (theta_1[j]-theta_0[j]) / 2);
                      }
#else
                      fb_mag = kb * tan( (theta_1[j]-theta_0[j]) / 2);
#endif

                      // each particle updates its force (no need for atomics)
                      int force_index = (i*NPC + j)*AMREX_SPACEDIM; 
                      fc_ptr[force_index  ] -= fb_mag*normal_sum_x;
                      fc_ptr[force_index+1] -= fb_mag*normal_sum_y;
                      
                      // These are the equal-and-opposite forces on the two ends of the line segment
                      int fi_1 = (i*NPC + j_mns)*AMREX_SPACEDIM; 
                      fc_ptr[fi_1  ] += 0.5*fb_mag*normal_sum_x;
                      fc_ptr[fi_1+1] += 0.5*fb_mag*normal_sum_y;
                      int fi_2 = (i*NPC + j_pls)*AMREX_SPACEDIM; 
                      fc_ptr[fi_2  ] += 0.5*fb_mag*normal_sum_x;
                      fc_ptr[fi_2+1] += 0.5*fb_mag*normal_sum_y;

                    }
                }        
            });

            // 
            /*********************************************************************
             * Calculate Fa (Area Force) -- SAME CELL ONLY                       * 
             ********************************************************************/
            amrex::ParallelFor(nrp,
              [pstruct,subdt,fc_ptr,nrp]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                ParticleType& p= pstruct[i];
                {
                    Real X_new[NPC] = {};
                    Real Y_new[NPC] = {};

                    Real cent_x(0.);
                    Real cent_y(0.);
 
                    // Here we iterate over all the vertices of this cell 
                    for (int j = 0; j < NPC; j++)
                    {
                        // Location of vertex j
                        Real x_0 = p.rdata(realData::pos_x0+j*AMREX_SPACEDIM);
                        Real y_0 = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);

                        // Velocity of vertex j
                        Real u_0 = p.rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM);
                        Real v_0 = p.rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM);

                        X_new[j] = x_0 + subdt * u_0;
                        Y_new[j] = y_0 + subdt * v_0;

                        cent_x += X_new[j];
                        cent_y += Y_new[j];
                    }

                    // Centroid of the polygon after it is evolve forward in time (temporarily)
                    cent_x /= NPC;
                    cent_y /= NPC;

                    Real circumference = 0.;

                    for (int j = 0; j < NPC; j++)
                    {
                        int j_mns = j;
                        int j_pls = (j == NPC -1) ? 0 : j+1;

                        Real x_mns = X_new[j_mns];
                        Real y_mns = Y_new[j_mns];

                        Real x_pls = X_new[j_pls];
                        Real y_pls = Y_new[j_pls];
   
                        Real side_length = std::sqrt( (x_pls-x_mns)*(x_pls-x_mns)+(y_pls-y_mns)*(y_pls-y_mns));
                        circumference += side_length;
                    }

                    Real a1 = areaOfPolygon(X_new,Y_new,NPC);

                    Real fa_mag = ka * (cell_a0 - a1) / (constant_r0*constant_r0);

                    for (int j = 0; j < NPC; j++)
                    {
                         int j_mns = j;
                         int j_pls = (j == NPC -1) ? 0 : j+1;

                         Real x_mns = X_new[j_mns];
                         Real y_mns = Y_new[j_mns];

                         Real x_pls = X_new[j_pls];
                         Real y_pls = Y_new[j_pls];

                         Real dist_x = x_pls-x_mns;
                         Real dist_y = y_pls-y_mns;
                         Real dist_mag = std::sqrt( dist_x*dist_x+dist_y*dist_y);
                         
                         Real side_length_frac = dist_mag / circumference;

                         RealVect normal(0.);
                         normal[0] =  dist_y / dist_mag;
                         normal[1] = -dist_x / dist_mag;

                         // We need to make sure the normal is *outward* facing relative to the centroid
                         // If the dot product of the (line from the centroid to the center of the line segment)
                         // with the (normal to the line segement) is negative then the normal is not outward-facing
                         // so we flip the sign
                         Real cent_to_seg_x = 0.5*(x_mns+x_pls) - cent_x;
                         Real cent_to_seg_y = 0.5*(y_mns+y_pls) - cent_y;
                         Real dot_product = cent_to_seg_x * normal[0] + cent_to_seg_y * normal[1];
                         if (dot_product < 0.) 
                         {
                            normal[0] *= -1.;
                            normal[1] *= -1.;
                         }

                         // Exert outward force on both ends of the line segment
                         int force_index_p = (i*NPC + j_mns)*AMREX_SPACEDIM; 
                         fc_ptr[force_index_p  ] += 0.5*fa_mag*normal[0]*side_length_frac;
                         fc_ptr[force_index_p+1] += 0.5*fa_mag*normal[1]*side_length_frac;

                         int force_index_m = (i*NPC + j_pls)*AMREX_SPACEDIM; 
                         fc_ptr[force_index_m  ] += 0.5*fa_mag*normal[0]*side_length_frac;
                         fc_ptr[force_index_m+1] += 0.5*fa_mag*normal[1]*side_length_frac;
                    }
                }
            });

            BL_PROFILE_VAR_STOP(calc_particle_collisions);
            BL_PROFILE_VAR("des::update_particle_velocity_and_position()", des_time_march);
            /********************************************************************
             * Update particle velocity and positions                           *
             *******************************************************************/

	    //Real sum(0.);

            amrex::ParallelFor(nrp,
              [pstruct,subdt,fc_ptr]
              AMREX_GPU_DEVICE (int i) noexcept
              //for (int i = 0; i < nrp; i++)
              {
                ParticleType& p = pstruct[i];
                {
                   // We define the location of the cell to be its centroid
                   Real cent_x(0.);
                   Real cent_y(0.);
            
                   /********************************************************************
                    * This block of code applies a force                               *
                    * to three particles with the greatest                             *
                    * local y value                                                    *
                    *******************************************************************/
#if 0
                   // Pick a cell to apply a constant force in input file
                   if (i == move_cell)
                   {
                      Real y_0_store[NPC] = {};

                      // Find the three largest y values in our cell
                      for (int j = 0; j < NPC; j++)
                      {
                         Real y_0_temp = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);
                         y_0_store[j] = y_0_temp;
                      }
                   
                      Real first, second, third;
                      Real first_index, second_index, third_index;
                      third = first = second = INT_MIN;

                      for(int i = 0; i <NPC; i++)
                      {
                         if (y_0_store[i] > first)
                         {
                            third = second;
                            second = first;
                            first = y_0_store[i];
                            first_index = i;
                         }
                         else if (y_0_store[i] > second)
                         {
                            third = second;
                            second = y_0_store[i];
                            second_index = i;
                         }
                         else if (y_0_store[i] > third)
                            third = y_0_store[i];
                            third_index = i;
                      }
                      //Real y_0_three[3] = {first, second, third}; // Store cell values
                      Real y_0_three_index[3] = {first_index, second_index, third_index}; // Store index of largest values

                      for (int k = 0; k < 3; k++)
                      {
                         int force_index_cell = (i*NPC + y_0_three_index[k])*AMREX_SPACEDIM; 
                         // Set to 0.0 and 0.0 to calculate from center
                         Real x_reference = 1.0;
                         Real y_reference = -2.0;

                         Real x_from_bottom = x_reference - p.pos(0);
                         Real y_from_bottom = y_reference - p.pos(1);
                         Real radius_from_reference = std::sqrt(x_from_bottom*x_from_bottom +y_from_bottom*y_from_bottom);

                         Real search_distance = 4.0;
                         if (radius_from_reference < search_distance)
                         {
                            // Linear constant Force
                            //fc_ptr[force_index_cell  ] += x_const_force * (4.0-radius_from_bottom);
                            //fc_ptr[force_index_cell+1] += y_const_force * (4.0-radius_from_bottom);
                            
                            // Log constant force
                            fc_ptr[force_index_cell  ] += x_const_force * (-1*std::log(radius_from_reference + 1)
                                  *0.5* radius_from_reference + search_distance);
                            fc_ptr[force_index_cell+1] += y_const_force * (-1*std::log(radius_from_reference + 1)
                                  *0.5* radius_from_reference + search_distance);
                         }
                      }
                   }
#endif

                   for (int j = 0; j < NPC; j++)
                   {
                       // These are offsets, not absolute positions 
                       Real x_0 = p.rdata(realData::pos_x0+j*AMREX_SPACEDIM);
                       Real y_0 = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);
                       Real u_0 = p.rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM);
                       Real v_0 = p.rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM);

                      /********************************************************************
                       * This block of code applies a force                               *
                       * to the entire cell. Comment this section                         *
                       * if applying forces to a select number of                         *
                       * vertices as above                                                *
                       *******************************************************************/
                       // Pick a cell to apply a constant force in input file
                       //if (i == move_cell)
                       //{
                       //   int force_index_cell = (i*NPC + j)*AMREX_SPACEDIM; 
                       //   Real radius_from_center = std::sqrt(p.pos(0)*p.pos(0) + p.pos(1)*p.pos(1));
                       //   if (radius_from_center < 1.0 && y_0 > 0)
                       //   {
                       //       fc_ptr[force_index_cell  ] += x_const_force * (1.0 - radius_from_center);
                       //       fc_ptr[force_index_cell+1] += y_const_force * (1.0 - radius_from_center);
                       //   }
                       //}
                       
#if 0
                       // Pick a cell to apply a constant force in input file
                       if (i == move_cell)
                       {
                          Real y_0_store[NPC] = {};
                          // Find the three largest y values in our cell
                          for (int j = 0; j < NPC; j++)
                          {
                             Real y_0_temp = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);
                             y_0_store[j] = y_0_temp;
                          }
                          Real first, second, third;
                          Real first_index, second_index, third_index;
                          third = first = second = INT_MIN;
                          for(int i = 0; i <NPC; i++)
                          {
                             if (y_0_store[i] > first)
                             {
                                third = second;
                                second = first;
                                first = y_0_store[i];
                                first_index = i;
                             }
                             else if (y_0_store[i] > second)
                             {
                                third = second;
                                second = y_0_store[i];
                                second_index = i;
                             }
                             else if (y_0_store[i] > third)
                                third = y_0_store[i];
                                third_index = i;
                          }
                          Real y_0_three[3] = {first, second, third}; // Store cell values
                          Real y_0_three_index[3] = {first_index, second_index, third_index}; // Store index of largest values
                          for (int k = 0; k < 3; k++)
                          {
                             int force_index_cell = (i*NPC + y_0_three_index[k])*AMREX_SPACEDIM;
                             Real radius_from_center = std::sqrt(p.pos(0)*p.pos(0) + p.pos(1)*p.pos(1));
                             //amrex::Print() <<  " cell : " << move_cell << std::endl;
                             //amrex::Print() <<  " y_0_three_index[k] : "<< y_0_three_index[k] << std::endl;
                             //amrex::Print() <<  " y_0_three[k] : "<< y_0_three[k] << std::endl;
                             //amrex::Print() <<  " radius_from_center: " << radius_from_center << std::endl;
                             if (radius_from_center < 3.0)
                             {
                               fc_ptr[force_index_cell  ] += x_const_force * (std::sqrt(3.0-radius_from_center));
                               fc_ptr[force_index_cell+1] += y_const_force * (std::sqrt(3.0-radius_from_center));
                             }
                          }
                       }
#endif

                       x_0 += subdt * u_0;
                       y_0 += subdt * v_0;
                       
                       int force_index = (i*NPC + j)*AMREX_SPACEDIM; 

                       // Up until now the forces have been scaled by constant_r0;
                       //    here we return them to physical units
                       u_0 += subdt * fc_ptr[force_index  ] * constant_r0;
                       v_0 += subdt * fc_ptr[force_index+1] * constant_r0;

                       //sum += subdt * fc_ptr[force_index+1];

                       // These are still offsets, not absolute positions 
                       p.rdata(realData::pos_x0+     j *AMREX_SPACEDIM) = x_0;
                       p.rdata(realData::pos_y0+     j *AMREX_SPACEDIM) = y_0;
                       p.rdata(realData::pos_x0+(NPC+j)*AMREX_SPACEDIM) = u_0;
                       p.rdata(realData::pos_y0+(NPC+j)*AMREX_SPACEDIM) = v_0;
                       
                       cent_x += x_0;
                       cent_y += y_0;
                   }

                   // This is the new average *offset*, not the actual centroid location in space
                   // Previous offset was 0 by construction
                   cent_x /= NPC;
                   cent_y /= NPC;

                   // This is real location in physical space
                   p.pos(0) += cent_x;
                   p.pos(1) += cent_y;

                   // Adjust the offsets so they are relative to the new value of the centroid
                   for (int j = 0; j < NPC; j++)
                   {
                       p.rdata(realData::pos_x0+     j *AMREX_SPACEDIM) -= cent_x;
                       p.rdata(realData::pos_y0+     j *AMREX_SPACEDIM) -= cent_y;
                   }

                }
              });
              // };

              //amrex::Print() << "SUM " << sum << std::endl; 

            BL_PROFILE_VAR_STOP(des_time_march);

            // ************************************************************************************************
            // Check for cell intersection -- does a vertex of one cell lie on the inside of another cell?    *
            // ************************************************************************************************
            amrex::ParallelFor(nrp,
                [pstruct,fc_ptr,nbor_data,nrp,subdt]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  Real X_poly[NPC+1] = {};
                  Real Y_poly[NPC+1] = {};

                  ParticleType& p1 = pstruct[i];
                  for (const auto& p2 : nbor_data.getNeighbors(i))
                  {
                      if (p1.idata(intData::ccell) !=  p2.idata(intData::ccell))
                      {
                        // We can fill X_poly and Y_poly before looping over j
                        for (int k = 0; k < NPC; k++)
                        {
                            // Location of vertices in polygon
                            X_poly[k] = p2.rdata(realData::pos_x0+k*AMREX_SPACEDIM) + p2.pos(0);
                            Y_poly[k] = p2.rdata(realData::pos_y0+k*AMREX_SPACEDIM) + p2.pos(1);
                        }

                        X_poly[NPC] = p2.rdata(realData::pos_x0+0*AMREX_SPACEDIM) + p2.pos(0);
                        Y_poly[NPC] = p2.rdata(realData::pos_y0+0*AMREX_SPACEDIM) + p2.pos(1);

                        // Here we iterate over all the vertices of cell p1 
                        for (int j = 0; j < NPC; j++)
                        {
                           // Location of vertex j of cell i
                           Real x_0 = p1.rdata(realData::pos_x0+j*AMREX_SPACEDIM) + p1.pos(0);
                           Real y_0 = p1.rdata(realData::pos_y0+j*AMREX_SPACEDIM) + p1.pos(1);
                             
                           bool check = cellPenetrationCheck(x_0, y_0, X_poly, Y_poly, NPC);

                           if ( check ){
                              amrex::Print() << "  " << std::endl;
                              amrex::Print() << " Vertex of cell 1 at " << RealVect(x_0,y_0) << std::endl;
                              amrex::Print() << " belonging to cell " 
                                             << p1.idata(intData::ccell) << " with centroid location " << RealVect(p1.pos(0),p1.pos(1))
                                             << std::endl;
                              amrex::Print() << "   is inside cell " << p2.idata(intData::ccell)
                                             << " with centroid location " << RealVect(p2.pos(0),p2.pos(1)) << std::endl;
                              amrex::Print() << "  " << std::endl;
                              Abort("\n!!! ");
                              }
                        }// vertices of p1
                      }// i != j
                  } // loop over p2
              }); // loop over p1

            Gpu::synchronize();

        } // particle tile

    } // end of computing forces, updating velocities and positions

    // ********************************************************************
    // Redistribute particles at the end of all substeps (note that the 
    // particle neighbour list needs to be reset when redistributing).
    // ********************************************************************
    if ((nstep+1) % nbor_interval == 0) {
        clearNeighbors();
        Redistribute(0, 0, 0, 1);
    }

    BL_PROFILE_REGION_STOP("ccm_dem::EvolveParticles()");
}
