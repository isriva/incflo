#include <ccm_pc.H>
#include <ccm_pc_K.H>

using namespace amrex;

void CCMParticleContainer::ComputeAndPrintDiagnostics (int lev,
                                                       int nstep,
                                                       Real time)
{
    BL_PROFILE_REGION_START("ccm_dem::ComputeDiagnostics()");
    BL_PROFILE("ccm_dem::ComputeDiagnostics()");

    /********************************************************************
     * Initialize area and length diagnostics
     *******************************************************************/
    Real area_min =  1.e30;
    Real area_max = -1.e30;
    Real area_tot = 0.;

    Real rad_min =  1.e30;
    Real rad_max = -1.e30;

    Real side_min =  1.e30;
    Real side_max = -1.e30;

    Real ratio_min =  1.e30;
    Real ratio_max = -1.e30;

    Real avg_x_centroid = 0.;
    Real avg_y_centroid = 0.;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (CCMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
            PairIndex index(pti.index(), pti.LocalTileIndex());

            const int nrp = GetParticles(lev)[index].numRealParticles();

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            // ************************************************************************************************
            // Calculate min area, max area and total area of all cells
            // ************************************************************************************************
         
            Real side;
            
            int bins[12] = {0};


            for (int i(0); i < nrp; ++i)
            {
                ParticleType& p= pstruct[i];
                
                Real X_new[NPC];
                Real Y_new[NPC];
 
                // We reset these for each polygon
                Real local_rad_min  =  1.e30; 
                Real local_rad_max  = -1.e30;
 
                // Initialize perimeter at 0
                Real perim = 0;

                // Here we iterate over all the vertices of this cell 
                for (int j = 0; j < NPC; j++)
                {
                    X_new[j] = p.rdata(realData::pos_x0+j*AMREX_SPACEDIM);
                    Y_new[j] = p.rdata(realData::pos_y0+j*AMREX_SPACEDIM);

                    // Remember that X and Y are relative to the centroid
                    Real rad = std::sqrt( X_new[j]*X_new[j] + Y_new[j]*Y_new[j]);

                    local_rad_min = std::min(local_rad_min, rad);
                    local_rad_max = std::max(local_rad_max, rad);

                    rad_min = std::min(rad_min, rad);
                    rad_max = std::max(rad_max, rad);


                    if (j > 0)
                    {
                        side = std::sqrt( (X_new[j]-X_new[j-1])*(X_new[j]-X_new[j-1]) + 
                                          (Y_new[j]-Y_new[j-1])*(Y_new[j]-Y_new[j-1]) );
                        side_min = std::min(side_min, side);
                        side_max = std::max(side_max, side);
                        perim = perim + side;
                    }
                
                } 

                // Now make sure we capture the last side
                side = std::sqrt( (X_new[0]-X_new[NPC-1])*(X_new[0]-X_new[NPC-1]) + 
                                  (Y_new[0]-Y_new[NPC-1])*(Y_new[0]-Y_new[NPC-1]) );
                side_min = std::min(side_min, side);
                side_max = std::max(side_max, side);
                perim = perim + side; 

                Real ratio = local_rad_max / local_rad_min;


# if 0
                if ( ratio <= 1.125 and ratio >= 1.000 ){
                   Bins[0] += 1;
                }
                else if ( ratio <= 1.250 and ratio > 1.125 ){
                   Bins[1] += 1;
                }
                else if ( ratio <= 1.375 and ratio > 1.250 ){
                   Bins[2] += 1;
                }
                else if ( ratio <= 1.500 and ratio > 1.375 ){
                   Bins[3] += 1;
                }
                else if ( ratio <= 1.625 and ratio > 1.500 ){
                   Bins[4] += 1;
                }
                else if ( ratio <= 1.750 and ratio > 1.625 ){
                   Bins[5] += 1;
                }
                else if ( ratio <= 1.875 and ratio > 1.750 ){
                   Bins[6] += 1;
                }
                else if ( ratio <= 2.000 and ratio > 1.875 ){
                   Bins[7] += 1;
                }
                else if ( ratio > 2.000 ){
                   Bins[8] += 1;
                }
# endif

                Real a1 = areaOfPolygon(X_new,Y_new,NPC);

                area_min = std::min(a1,area_min);
                area_max = std::max(a1,area_max);
                area_tot += a1;

                ratio_min = std::min(ratio,ratio_min);
                ratio_max = std::max(ratio,ratio_max);

                avg_x_centroid += p.pos(0);
                avg_y_centroid += p.pos(1);

                //calculate asphericity for each cell
                Real asph = ((perim*perim)/(4*M_PI*a1));
                //amrex::Print() << " asphericity: " << asph << std::endl;  

                if (asph >= 1.0 && asph < 1.05)
                    bins[0] += 1;
                else if (asph >= 1.05 && asph < 1.1)
                    bins[1] += 1;
                else if (asph >= 1.1 && asph < 1.15)
                    bins[2] += 1;
                else if (asph >= 1.15 && asph < 1.2)
                    bins[3] += 1;
                else if (asph >= 1.2 && asph < 1.25)
                    bins[4] += 1;
                else if (asph >= 1.25 && asph < 1.3)
                    bins[5] += 1;
                else if (asph >= 1.3 && asph < 1.35)
                    bins[6] += 1;
                else if (asph >= 1.35 && asph < 1.4)
                    bins[7] += 1;
                else if (asph >= 1.4 && asph < 1.45)
                    bins[8] += 1;
                else if (asph >= 1.45 && asph < 1.5)
                    bins[9] += 1;
                else if (asph >= 1.5 && asph < 1.55)
                    bins[10] += 1;
                else if (asph >= 1.55 && asph < 1.6)
                    bins[11] += 1;
            }
            int binsTotal = 0;
            for (int l = 0; l < 12; l++)
            binsTotal += bins[l];
            amrex::Print() << " 1.0-1.05: " << bins[0] << std::endl;                              
            amrex::Print() << " 1.05-1.1: " << bins[1] << std::endl;
            amrex::Print() << " 1.1-1.15: " << bins[2] << std::endl;
            amrex::Print() << " 1.15-1.2: " << bins[3] << std::endl;
            amrex::Print() << " 1.2-1.25: " << bins[4] << std::endl;
            amrex::Print() << " 1.25-1.3: " << bins[5] << std::endl;
            amrex::Print() << " 1.3-1.35: " << bins[6] << std::endl;
            amrex::Print() << " 1.35-1.4: " << bins[7] << std::endl;
            amrex::Print() << " 1.4-1.45: " << bins[8] << std::endl;
            amrex::Print() << " 1.45-1.5: " << bins[9] << std::endl;
            amrex::Print() << " 1.5-1.55: " << bins[10] << std::endl;
            amrex::Print() << " 1.55-1.6: " << bins[11] << std::endl;
            amrex::Print() << " binsTotal: " << binsTotal << std::endl;

# if 0
            int binsTotal = 0;
            for ( int i = 0; i < 9; i++ ){
               binsTotal += Bins[i];
            }
            amrex::Print() << " Bins[0]: " << Bins[0] << std::endl; 
            amrex::Print() << " Bins[1]: " << Bins[1] << std::endl; 
            amrex::Print() << " Bins[2]: " << Bins[2] << std::endl; 
            amrex::Print() << " Bins[3]: " << Bins[3] << std::endl; 
            amrex::Print() << " Bins[4]: " << Bins[4] << std::endl; 
            amrex::Print() << " Bins[5]: " << Bins[5] << std::endl; 
            amrex::Print() << " Bins[6]: " << Bins[6] << std::endl; 
            amrex::Print() << " Bins[7]: " << Bins[7] << std::endl; 
            amrex::Print() << " Bins[8]: " << Bins[8] << std::endl; 
            amrex::Print() << " binsTotal: " << binsTotal << std::endl; 
# endif
    }


    /********************************************************************
     * Finish the area diagnostics by parallel reduce, scaling and printing
     *******************************************************************/
    ParallelDescriptor::ReduceRealSum(area_tot); 
    ParallelDescriptor::ReduceRealMin(area_min); 
    ParallelDescriptor::ReduceRealMax(area_max); 

    ParallelDescriptor::ReduceRealSum(avg_x_centroid); 
    ParallelDescriptor::ReduceRealSum(avg_y_centroid); 
    long total_num = TotalNumberOfParticles();
    avg_x_centroid /= static_cast<Real>(total_num);
    avg_y_centroid /= static_cast<Real>(total_num);

    ParallelDescriptor::ReduceRealMin(rad_min); 
    ParallelDescriptor::ReduceRealMax(rad_max); 

    ParallelDescriptor::ReduceRealMin(side_min); 
    ParallelDescriptor::ReduceRealMax(side_max); 

    ParallelDescriptor::ReduceRealMin(ratio_min); 
    ParallelDescriptor::ReduceRealMax(ratio_max); 

    amrex::Print() << "Min / max / average fraction of cell areas: " 
                   << area_min/cell_a0 << " " << area_max/cell_a0 << " " << area_tot/orig_domain_area << std::endl;
    amrex::Print() << "Min rad  / Max rad  relative to cell_r0   : "
                   << rad_min/cell_r0 << " " << rad_max/cell_r0 << std::endl;
    amrex::Print() << "Min side / Max side relative to cell_x0   : "
                   << side_min/cell_x0 << " " << side_max/cell_x0 << std::endl;
    amrex::Print() << "Min ratio / Max ratio over all cells      : "
                   << ratio_min << " " << ratio_max << std::endl;

    cell_rmax = rad_max;

    std::ofstream runlog;
    std::string filename = "out_area"; 
    runlog.open(filename.c_str(),std::ios::out|std::ios::app);
    if (!runlog.good())
        amrex::FileOpenFailed(filename);
   
    // Only write every 10 steps
    if (nstep%10 == 0)
       runlog << std::setw(8) << time << " " << area_tot/orig_domain_area << std::endl;

    std::ofstream runlog2;
    std::string filename2 = "out_centroid"; 
    runlog2.open(filename2.c_str(),std::ios::out|std::ios::app);
    if (!runlog2.good())
        amrex::FileOpenFailed(filename2);
   
    // Only write every 10 steps
    if (nstep%10 == 0)
       runlog2 << std::setw(8) << time << " " << avg_x_centroid << " " << avg_y_centroid << std::endl;

    BL_PROFILE_REGION_STOP("ccm_dem::ComputeDiagnostics()");
}
