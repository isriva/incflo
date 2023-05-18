#ifndef _TILE_NO_Z_H_
#define _TILE_NO_Z_H_

#include <AMReX.H>

/**
 * Function returns explicit tile size that prevents tiling in z
 */

AMREX_FORCE_INLINE
amrex::IntVect TileNoZ()
{
    if (amrex::TilingIfNotGPU()) {
#if (AMREX_SPACEDIM==3)
        return amrex::IntVect{amrex::FabArrayBase::mfiter_tile_size[0],
                              amrex::FabArrayBase::mfiter_tile_size[1],
                              1024000};
#elif (AMREX_SPACEDIM==2)
        return amrex::IntVect{amrex::FabArrayBase::mfiter_tile_size[0],
                              1024000};
#endif
    }
    else {
        return amrex::IntVect::TheZeroVector();
    }
}

#endif
