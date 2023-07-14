#include <incflo.H>

using namespace amrex;

void incflo::prob_init_fluid (int lev)
{
    auto& ld = *m_leveldata[lev];
    Box const& domain = geom[lev].Domain();
    auto const& dx = geom[lev].CellSizeArray();
    auto const& problo = geom[lev].ProbLoArray();
    auto const& probhi = geom[lev].ProbHiArray();

    ld.p_nd.setVal(0.0);
    ld.gp.setVal(0.0);
    ld.vel_eta.setVal(1.0);
    ld.vel_eta1.setVal(1.0);
    ld.vel_eta2.setVal(1.0);

    ld.density.setVal(m_ro_0);
    ld.density_o.setVal(m_ro_0);
    ld.density0.setVal(0.0);

    AMREX_D_TERM(ld.velocity.setVal(m_ic_u, 0, 1);,
                 ld.velocity.setVal(m_ic_v, 1, 1);,
                 ld.velocity.setVal(m_ic_w, 2, 1););

    if (m_ntrac > 0) ld.tracer.setVal(0.0);

    for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
//        const Box& vbx = mfi.growntilebox(1);
        const Box& gbx = mfi.fabbox();
        const Box& nbx = mfi.tilebox(IntVect::TheNodeVector()); // nodal box
        if (0 == m_probtype || 114 == m_probtype )
        { }
        else if (1 == m_probtype)
        {
            init_taylor_green(vbx, gbx,
                              ld.velocity.array(mfi),
                              ld.density.array(mfi),
                              ld.tracer.array(mfi),
                              domain, dx, problo, probhi);
        }
        else if (2 == m_probtype)
        {
            init_taylor_vortex(vbx, gbx,
                               ld.velocity.array(mfi),
                               ld.density.array(mfi),
                               ld.tracer.array(mfi),
                               domain, dx, problo, probhi);
        }
        else if (3 == m_probtype)
        {
            init_taylor_green3d(vbx, gbx,
                                ld.velocity.array(mfi),
                                ld.density.array(mfi),
                                ld.tracer.array(mfi),
                                domain, dx, problo, probhi);
        }
        else if (4 == m_probtype)
        {
            init_couette(vbx, gbx,
                         ld.velocity.array(mfi),
                         ld.density.array(mfi),
                         ld.tracer.array(mfi),
                         domain, dx, problo, probhi);
        }
        else if (5 == m_probtype)
        {
            init_rayleigh_taylor(vbx, gbx,
                                 ld.velocity.array(mfi),
                                 ld.density.array(mfi),
                                 ld.tracer.array(mfi),
                                 domain, dx, problo, probhi);
        }
        else if (6 == m_probtype)
        {
            init_channel_slant(vbx, gbx,
                               ld.velocity.array(mfi),
                               ld.density.array(mfi),
                               ld.tracer.array(mfi),
                               domain, dx, problo, probhi);
        }
        else if (11 == m_probtype)
        {
            init_tuscan(vbx, gbx,
                        ld.velocity.array(mfi),
                        ld.density.array(mfi),
                        ld.tracer.array(mfi),
                        domain, dx, problo, probhi);
        }
        else if (111 == m_probtype || 112 == m_probtype || 113 == m_probtype)
        {
            init_boussinesq_bubble(vbx, gbx,
                                   ld.velocity.array(mfi),
                                   ld.density.array(mfi),
                                   ld.tracer.array(mfi),
                                   domain, dx, problo, probhi);
        }
        else if (12 == m_probtype)
        {
            init_periodic_tracer(vbx, gbx,
                                 ld.velocity.array(mfi),
                                 ld.density.array(mfi),
                                 ld.tracer.array(mfi),
                                 domain, dx, problo, probhi);
        }
        else if (13 == m_probtype)
        {
            init_flow_in_box(vbx, gbx,
                             ld.velocity.array(mfi),
                             ld.density.array(mfi),
                             ld.tracer.array(mfi),
                             domain, dx, problo, probhi);
        }
        else if (14 == m_probtype)
        {
            init_circ_traceradvect(vbx, gbx,
                                   ld.velocity.array(mfi),
                                   ld.density.array(mfi),
                                   ld.tracer.array(mfi),
                                   domain, dx, problo, probhi);
        }
        else if (15 == m_probtype)
        {
            init_gaussian_traceradvect(vbx, gbx,
                                       ld.velocity.array(mfi),
                                       ld.density.array(mfi),
                                       ld.tracer.array(mfi),
                                       domain, dx, problo, probhi);
        }
        else if (66 == m_probtype)
        {
            init_vortex_in_sphere(vbx, gbx,
                                  ld.velocity.array(mfi),
                                  ld.density.array(mfi),
                                  ld.tracer.array(mfi),
                                  domain, dx, problo, probhi);
        }
        else if (21 == m_probtype || 22 == m_probtype || 23 == m_probtype)
        {
            init_double_shear_layer(vbx, gbx,
                                    ld.velocity.array(mfi),
                                    ld.density.array(mfi),
                                    ld.tracer.array(mfi),
                                    domain, dx, problo, probhi);
        }
        else if (31  == m_probtype || 32  == m_probtype || 33  == m_probtype ||
                 311 == m_probtype || 322 == m_probtype || 333 == m_probtype ||
                 41  == m_probtype)
        {
            init_plane_poiseuille(vbx, gbx,
                                  ld.velocity.array(mfi),
                                  ld.density.array(mfi),
                                  ld.tracer.array(mfi),
                                  domain, dx, problo, probhi);
        }
        else if (51 == m_probtype)
        {
            init_rayleigh_taylor_vof(vbx, gbx,
                                     ld.velocity.array(mfi),
                                     ld.density.array(mfi),
                                     ld.tracer.array(mfi),
                                     domain, dx, problo, probhi);
        }
        else if (52 == m_probtype)
        {
            inclined_channel(vbx, gbx, nbx, 
                            ld.velocity.array(mfi),
                            ld.density.array(mfi),
                            ld.tracer.array(mfi),
                            ld.gp.array(mfi),
                            domain, dx, problo, probhi);
        }
        else if (521 == m_probtype)
        {
            inclined_channel_granular(vbx, gbx, nbx, 
                                      ld.velocity.array(mfi),
                                      ld.density.array(mfi),
                                      ld.tracer.array(mfi),
                                      ld.p_nd.array(mfi),
                                      ld.p0.array(mfi),
                                      ld.p_visc.array(mfi),
                                      ld.gp.array(mfi),
                                      ld.gp0.array(mfi),
                                      ld.density0.array(mfi),
                                      domain, dx, problo, probhi,0);
        }
        else if (522 == m_probtype)
        {
            inclined_channel_granular(vbx, gbx, nbx, 
                                      ld.velocity.array(mfi),
                                      ld.density.array(mfi),
                                      ld.tracer.array(mfi),
                                      ld.p_nd.array(mfi),
                                      ld.p0.array(mfi),
                                      ld.p_visc.array(mfi),
                                      ld.gp.array(mfi),
                                      ld.gp0.array(mfi),
                                      ld.density0.array(mfi),
                                      domain, dx, problo, probhi,1);
        }
        else if (53 == m_probtype)
        {
            column_collapse (vbx, gbx, nbx, 
                            ld.velocity.array(mfi),
                            ld.density.array(mfi),
                            ld.tracer.array(mfi),
                            ld.gp.array(mfi),
                            ld.p0.array(mfi),
                            domain, dx, problo, probhi);
        }
        else if (531 == m_probtype)
        {
            column_collapse_granular (vbx, gbx, nbx, 
                                      ld.velocity.array(mfi),
                                      ld.density.array(mfi),
                                      ld.tracer.array(mfi),
                                      ld.p_nd.array(mfi),
                                      ld.p0.array(mfi),
                                      ld.p_visc.array(mfi),
                                      ld.gp.array(mfi),
                                      ld.gp0.array(mfi),
                                      ld.density0.array(mfi),
                                      domain, dx, problo, probhi);
        }
#if 0
        else if (500 == m_probtype)
        {
            init_heated_ground(vbx, gbx,
                               ld.velocity.array(mfi),
                               ld.density.array(mfi),
                               ld.tracer.array(mfi),
                               domain, dx, problo, probhi);
        }
#endif
        else
        {
            amrex::Abort("prob_init_fluid: unknown m_probtype");
        };
    }
}

void incflo::init_taylor_green (Box const& vbx, Box const& /*gbx*/,
                                Array4<Real> const& vel,
                                Array4<Real> const& /*density*/,
                                Array4<Real> const& /*tracer*/,
                                Box const& /*domain*/,
                                GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = Real(i+0.5)*dx[0];
        Real y = Real(j+0.5)*dx[1];
        constexpr Real twopi = Real(2.0)*Real(3.1415926535897932);
        vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y);
        vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y);
#if (AMREX_SPACEDIM == 3)
        vel(i,j,k,2) = Real(0.0);
#endif
    });
}

void incflo::init_taylor_green3d (Box const& vbx, Box const& /*gbx*/,
                                  Array4<Real> const& vel,
                                  Array4<Real> const& /*density*/,
                                  Array4<Real> const& /*tracer*/,
                                  Box const& /*domain*/,
                                  GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                  GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                  GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = Real(i+0.5)*dx[0];
        Real y = Real(j+0.5)*dx[1];
        Real z = Real(k+0.5)*dx[2];
        constexpr Real twopi = Real(2.0)*Real(3.1415926535897932);
        AMREX_D_TERM(vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y) * std::cos(twopi*z);,
                     vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y) * std::cos(twopi*z);,
                     vel(i,j,k,2) = 0.0;);
    });
}

void incflo::init_taylor_vortex (Box const& vbx, Box const& /*gbx*/,
                                 Array4<Real> const& vel,
                                 Array4<Real> const& /*density*/,
                                 Array4<Real> const& /*tracer*/,
                                 Box const& /*domain*/,
                                 GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                 GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                 GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = Real(i+0.5)*dx[0];
        Real y = Real(j+0.5)*dx[1];
        constexpr Real pi = Real(3.1415926535897932);
        constexpr Real u0 = Real(1.0);
        constexpr Real v0 = Real(1.0);
        vel(i,j,k,0) =  u0 - std::cos(pi*x) * std::sin(pi*y);
        vel(i,j,k,1) =  v0 + std::sin(pi*x) * std::cos(pi*y);
#if (AMREX_SPACEDIM == 3)
        vel(i,j,k,2) = Real(0.0);
#endif
    });
}


void incflo::init_vortex_in_sphere (Box const& vbx, Box const& /*gbx*/,
                               Array4<Real> const& vel,
                               Array4<Real> const& /*density*/,
                               Array4<Real> const& /*tracer*/,
                               Box const& /*domain*/,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               GpuArray<Real, AMREX_SPACEDIM> const& problo,
                               GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = problo[0] + Real(i+0.5)*dx[0];
        Real y = problo[1] + Real(j+0.5)*dx[1];
        Real deltax = x;
        Real deltay = y;
        Real d_sq = deltax*deltax + deltay*deltay;
        Real r_sq = Real(0.003)*Real(0.003);
        Real u_vort = -Real(0.2)*deltay/r_sq * std::exp(-d_sq/r_sq/Real(2.0));
        Real v_vort =  Real(0.2)*deltax/r_sq * std::exp(-d_sq/r_sq/Real(2.0));
        vel(i,j,k,0) =  u_vort;
        vel(i,j,k,1) =  v_vort;
#if (AMREX_SPACEDIM == 3)
        vel(i,j,k,2) = Real(0.0);
#endif
    });
}

void incflo::init_flow_in_box (Box const& vbx, Box const& /*gbx*/,
                               Array4<Real> const& vel,
                               Array4<Real> const& /*density*/,
                               Array4<Real> const& tracer,
                               Box const& /*domain*/,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                               GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    ParmParse pp("box");

    Vector<Real> boxLo(AMREX_SPACEDIM), boxHi(AMREX_SPACEDIM);
    Real offset = Real(1.0e-15);

    for(int i = 0; i < AMREX_SPACEDIM; i++)
    {
        boxLo[i] = geom[0].ProbLo(i);
        boxHi[i] = geom[0].ProbHi(i);
    }

    pp.queryarr("Lo", boxLo, 0, AMREX_SPACEDIM);
    pp.queryarr("Hi", boxHi, 0, AMREX_SPACEDIM);

#if (AMREX_SPACEDIM == 3)
    int periodic_dir;
    pp.get("periodic_dir", periodic_dir);
#endif

    pp.query("offset", offset);

    Real xlo = boxLo[0] + offset;
    Real xhi = boxHi[0] - offset;

    Real ylo = boxLo[1] + offset;
    Real yhi = boxHi[1] - offset;

#if (AMREX_SPACEDIM == 3)
    Real zlo = boxLo[2] + offset;
    Real zhi = boxHi[2] - offset;
#endif

#if (AMREX_SPACEDIM == 3)
    if (periodic_dir == 0)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = Real(j+0.5)*dx[1]*(yhi-ylo) + ylo;
            Real z = Real(k+0.5)*dx[2]*(zhi-zlo) + zlo;
            constexpr Real pi = Real(3.1415926535897932);
            vel(i,j,k,1) =  std::sin(pi*y) * std::cos(pi*z);
            vel(i,j,k,2) = -std::cos(pi*y) * std::sin(pi*z);
            vel(i,j,k,0) = Real(1.0);
            if (y < 0.5)
                tracer(i,j,k) = Real(0.0);
            else
                tracer(i,j,k) = Real(1.0);
        });
    } else if (periodic_dir == 1)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5)*dx[0]*(xhi-xlo) + xlo;
            Real z = Real(k+0.5)*dx[2]*(zhi-zlo) + zlo;
            constexpr Real pi = Real(3.1415926535897932);
            vel(i,j,k,2) =  std::sin(pi*z) * std::cos(pi*x);
            vel(i,j,k,0) = -std::cos(pi*z) * std::sin(pi*x);
            vel(i,j,k,1) = Real(1.0);
            if (z < Real(0.5))
                tracer(i,j,k) = Real(0.0);
            else
                tracer(i,j,k) = Real(1.0);
        });
    } else if (periodic_dir == 2)
#endif
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5)*dx[0]*(xhi-xlo) + xlo;
            Real y = Real(j+0.5)*dx[1]*(yhi-ylo) + ylo;
            constexpr Real pi = Real(3.1415926535897932);
            vel(i,j,k,0) =  std::sin(pi*x) * std::cos(pi*y);
            vel(i,j,k,1) = -std::cos(pi*x) * std::sin(pi*y);
#if (AMREX_SPACEDIM == 3)
            vel(i,j,k,2) = Real(1.0);
#endif
            if (x < Real(0.5))
                tracer(i,j,k) = Real(0.0);
            else
                tracer(i,j,k) = Real(1.0);
        });
#if (AMREX_SPACEDIM == 3)
    } else {
       amrex::Error("flow_in_box assumes a periodic direction");
#endif
    }
}

void incflo::init_circ_traceradvect (Box const& vbx, Box const& /*gbx*/,
                                     Array4<Real> const& vel,
                                     Array4<Real> const& density,
                                     Array4<Real> const& tracer,
                                     Box const& /*domain*/,
                                     GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                     GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                     GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];

        vel(i,j,k,0) = 1.;
        vel(i,j,k,1) = 0.5;

        density(i,j,k) = 1.;

        Real sum = 0.;
        for (int jj=0; jj<10; ++jj) {
            Real yy = (j + (jj+0.5)/10.) * dx[1];
            for (int ii=0; ii<10; ++ii) {
                Real xx = (i + (ii+0.5)/10.) * dx[0];

                Real r = std::sqrt( (xx-0.5)*(xx-0.5) + (yy-0.5)*(yy-0.5) );

                if (r < 0.1) {
                    sum += 1.;
                } else if (r == 0.1) {
                    sum += 0.5;
                }

            }
        }

        tracer(i,j,k) = sum / 100.;

    });

#elif (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        Real z = (k+0.5)*dx[2];

        vel(i,j,k,0) = 1.;
        vel(i,j,k,1) = 0.5;
        vel(i,j,k,2) = 0.25;

        density(i,j,k) = 1.;

        Real sum = 0.;
        for (int kk=0; kk<10; ++kk) {
            Real zz = (k + (kk+0.5)/10.) * dx[2];
            for (int jj=0; jj<10; ++jj) {
                Real yy = (j + (jj+0.5)/10.) * dx[1];
                for (int ii=0; ii<10; ++ii) {
                    Real xx = (i + (ii+0.5)/10.) * dx[0];

                    Real r = std::sqrt( (xx-0.5)*(xx-0.5) + (yy-0.5)*(yy-0.5) + (zz-0.5)*(zz-0.5) );

                    if (r < 0.1) {
                        sum += 1.;
                    } else if (r == 0.1) {
                        sum += 0.5;
                    }

                }
            }
        }

        tracer(i,j,k) = sum / 1000.;

    });
#endif

}void incflo::init_gaussian_traceradvect (Box const& vbx, Box const& /*gbx*/,
                                          Array4<Real> const& vel,
                                          Array4<Real> const& density,
                                          Array4<Real> const& tracer,
                                          Box const& /*domain*/,
                                          GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                          GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                          GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];


        vel(i,j,k,0) = 1.;
        vel(i,j,k,1) = 1.;

        density(i,j,k) = 1.;

#if (AMREX_SPACEDIM == 2)
        Real r = std::sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
#elif (AMREX_SPACEDIM == 3)
        Real z = (k+0.5)*dx[2];
        vel(i,j,k,2) = 1.;

        Real r = std::sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5) );
#endif

        tracer(i,j,k) = exp(-300.*r*r);
    });

}

void incflo::init_couette (Box const& vbx, Box const& /*gbx*/,
                           Array4<Real> const& vel,
                           Array4<Real> const& /*density*/,
                           Array4<Real> const& /*tracer*/,
                           Box const& domain,
                           GpuArray<Real, AMREX_SPACEDIM> const& /*dx*/,
                           GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                           GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    Real num_cells_y = static_cast<Real>(domain.length(1));
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = Real(j+0.5) / num_cells_y;
        AMREX_D_TERM(vel(i,j,k,0) *= (y-Real(0.5));,
                     vel(i,j,k,1) = Real(0.0);,
                     vel(i,j,k,2) = Real(0.0););
    });
}

void incflo::init_channel_slant (Box const& vbx, Box const& /*gbx*/,
                                 Array4<Real> const& /*vel*/,
                                 Array4<Real> const& density,
                                 Array4<Real> const& tracer,
                                 Box const& domain,
                                 GpuArray<Real, AMREX_SPACEDIM> const& /*dx*/,
                                 GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                 GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    const auto dhi = amrex::ubound(domain);
    int direction  = -1;

    // Get cylinder information from inputs file.                               *
    ParmParse pp("cylinder");
    pp.get("direction",  direction);

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (density(i,j,k)>0) {

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n)
                tracer(i,j,k,n) = Real(0.0);

            if (direction == 0) {
                if (nt > 0 && i <= dhi.x/8)   tracer(i,j,k,0) = Real(1.0);
                if (nt > 1 && i <= dhi.x/2)   tracer(i,j,k,1) = Real(2.0);
                if (nt > 2 && i <= dhi.x*3/4) tracer(i,j,k,2) = Real(3.0);
            } else if (direction == 1) {
                if (nt > 0 && j <= dhi.y/8)   tracer(i,j,k,0) = Real(1.0);
                if (nt > 1 && j <= dhi.y/2)   tracer(i,j,k,1) = Real(2.0);
                if (nt > 2 && j <= dhi.y*3/4) tracer(i,j,k,2) = Real(3.0);
#if (AMREX_SPACEDIM == 3)
            } else {
                if (nt > 0 && k <= dhi.z/8)   tracer(i,j,k,0) = Real(1.0);
                if (nt > 1 && k <= dhi.z/2)   tracer(i,j,k,1) = Real(2.0);
                if (nt > 2 && k <= dhi.z*3/4) tracer(i,j,k,2) = Real(3.0);
#endif
            }
        }
    });
}

void incflo::init_rayleigh_taylor (Box const& vbx, Box const& /*gbx*/,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& density,
                                   Array4<Real> const& tracer,
                                   Box const& /*domain*/,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    static constexpr Real pi    = Real(3.1415926535897932);
    static constexpr Real rho_1 = Real(0.5);
    static constexpr Real rho_2 = Real(2.0);
    static constexpr Real tra_1 = Real(0.0);
    static constexpr Real tra_2 = Real(1.0);

    static constexpr Real width = Real(0.005);

    const Real splitx = Real(0.5)*(problo[0] + probhi[0]);
    const Real L_x    = probhi[0] - problo[0];

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,0) = Real(0.0);
        vel(i,j,k,1) = Real(0.0);

        Real x = problo[0] + Real(i+0.5)*dx[0];
        Real y = problo[1] + Real(j+0.5)*dx[1];

        const Real r2d = amrex::min(std::abs(x-splitx), Real(0.5)*L_x);
        const Real pertheight = Real(0.5) - Real(0.01)*std::cos(Real(2.0)*pi*r2d/L_x);

        density(i,j,k) = rho_1 + (Real(0.5)*(rho_2-rho_1))*(Real(1.0)+std::tanh((y-pertheight)/width));
        tracer(i,j,k)  = tra_1 + (Real(0.5)*(tra_2-tra_1))*(Real(1.0)+std::tanh((y-pertheight)/width));
    });

#elif (AMREX_SPACEDIM == 3)

    const Real splity = Real(0.5)*(problo[1] + probhi[1]);
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,0) = Real(0.0);
        vel(i,j,k,1) = Real(0.0);
        vel(i,j,k,2) = Real(0.0);

        Real x = problo[0] + Real(i+0.5)*dx[0];
        Real y = problo[1] + Real(j+0.5)*dx[1];
        Real z = problo[2] + Real(k+0.5)*dx[2];

        const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), Real(0.5)*L_x);
        const Real pertheight = Real(0.5) - Real(0.01)*std::cos(Real(2.0)*pi*r2d/L_x);

        density(i,j,k) = rho_1 + (Real(0.5)*(rho_2-rho_1))*(Real(1.0)+std::tanh((z-pertheight)/width));
        tracer(i,j,k)  = tra_1 + (Real(0.5)*(tra_2-tra_1))*(Real(1.0)+std::tanh((z-pertheight)/width));
    });
#endif
}

void incflo::init_rayleigh_taylor_vof (Box const& vbx, Box const& /*gbx*/,
                                       Array4<Real> const& vel,
                                       Array4<Real> const& density,
                                       Array4<Real> const& tracer,
                                       Box const& /*domain*/,
                                       GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                       GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                       GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    static constexpr Real pi   = 3.1415926535897932;
    Real rho_1 = m_fluid_vof[0].rho;
    Real rho_2 = m_fluid_vof[1].rho;
    amrex::Print() << "rho1 and rho2 during RT setup: " << rho_1 << " " << rho_2 << std::endl;
    static constexpr Real tra_1 = 1.0;
    static constexpr Real tra_2 = 0.0;

    static constexpr Real width = 0.005;

    const Real splitx = 0.5*(problo[0] + probhi[0]);
    const Real L_x    = probhi[0] - problo[0];

#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;

        Real x = problo[0] + (i+Real(0.5))*dx[0];
        Real y = problo[1] + (j+Real(0.5))*dx[1];

        const Real r2d = amrex::min(std::abs(x-splitx), Real(0.5)*L_x);
        const Real pertheight = 0.5 - 0.01*std::cos(2.0*pi*r2d/L_x);

        Real conc = tra_1 + ((tra_2-tra_1)/2.0)*(1.0+std::tanh((y-pertheight)/width));
        density(i,j,k) = (rho_1*rho_2)/(conc*rho_2 + (1.0-conc)*rho_1);
//        density(i,j,k) = rho_1 + ((rho_2-rho_1)/2.0)*(1.0+std::tanh((y-pertheight)/width));
        tracer(i,j,k)  = conc;
    });

#elif (AMREX_SPACEDIM == 3)

    const Real splity = Real(0.5)*(problo[1] + probhi[1]);
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;

        Real x = problo[0] + (i+0.5)*dx[0];
        Real y = problo[1] + (j+0.5)*dx[1];
        Real z = problo[2] + (k+0.5)*dx[2];

        const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), Real(0.5)*L_x);
        const Real pertheight = Real(0.5) - 0.01*std::cos(Real(2.0)*pi*r2d/L_x);

        Real conc = tra_1 + ((tra_2-tra_1)/2.0)*(1.0+std::tanh((z-pertheight)/width));
        density(i,j,k) = (rho_1*rho_2)/(conc*rho_2 + (1.0-conc)*rho_1);
        tracer(i,j,k)  = conc;

        // density(i,j,k) = rho_1 + ((rho_2-rho_1)/2.0)*(1.0+std::tanh((z-pertheight)/width));
        // tracer(i,j,k)  = tra_1 + ((tra_2-tra_1)/2.0)*(1.0+std::tanh((z-pertheight)/width));
    });
#endif
}

void incflo::column_collapse  (Box const& vbx, Box const& /*gbx*/, Box const& nbx,
                               Array4<Real> const& vel,
                               Array4<Real> const& density,
                               Array4<Real> const& tracer,
                               Array4<Real> const& gp,
                               Array4<Real> const& p0,
                               Box const& /*domain*/,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               GpuArray<Real, AMREX_SPACEDIM> const& problo,
                               GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real rho_1, rho_2, rho;
    if (m_do_vof) {
        rho_1 = m_fluid_vof[0].rho;
        rho_2 = m_fluid_vof[1].rho;
        amrex::Print() << "rho1 and rho2 during column collapse setup: " << rho_1 << " " << rho_2 << std::endl;
    }
    else {
        amrex::Abort("need vof for this setup");
    }
    if ((m_gran_lim[0] < Real(0.0)) or (m_gran_lim[0] < Real(0.0))) amrex::Abort("provide m_gran_lim for this problem");

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
        Real y = problo[1] + (j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + (k+0.5)*dx[2];

        if ((x < m_gran_lim[0]) and (z < m_gran_lim[1])) {
            density(i,j,k) = rho_2;
            tracer(i,j,k) = 1.0;
        }
        else {
            density(i,j,k) = rho_1;
            tracer(i,j,k) = 0.0;
        }
        gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
        gp(i,j,k,1) = 0.0;
        gp(i,j,k,2) = m_gravity[2] * density(i,j,k);

        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
        if ((x < m_gran_lim[0]) and (y < m_gran_lim[1])) {
            density(i,j,k) = rho_2;
            tracer(i,j,k) = 1.0;
        }
        else {
            density(i,j,k) = rho_1;
            tracer(i,j,k) = 0.0;
        }
        gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
        gp(i,j,k,1) = m_gravity[1] * density(i,j,k);

        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
#endif

    });
}
void incflo::column_collapse_granular  (Box const& vbx, Box const& /*gbx*/, Box const& nbx,
                                        Array4<Real> const& vel,
                                        Array4<Real> const& density,
                                        Array4<Real> const& tracer,
                                        Array4<Real> const& p_nd,
                                        Array4<Real> const& p0,
                                        Array4<Real> const& p_visc,
                                        Array4<Real> const& gp,
                                        Array4<Real> const& gp0,
                                        Array4<Real> const& density0,
                                        Box const& /*domain*/,
                                        GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                        GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                        GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    if (!m_do_vof) amrex::Abort("need vof for this setup");
    if ((m_gran_lim[0] < Real(0.0)) or (m_gran_lim[1] < Real(0.0))) amrex::Abort("provide m_gran_lim for this problem");
    
    Real rho_1, rho_2, rho;
    rho_1 = m_fluid_vof[0].rho;
    rho_2 = m_fluid_vof[1].rho;
    amrex::Print() << "rho1 and rho2 during column collapse setup: " << rho_1 << " " << rho_2 << std::endl;

#if (AMREX_SPACEDIM == 3)
    const Real H = probhi[2] - problo[2];
    amrex::Print() << "granular region: " << problo[0] << " < x < " << m_gran_lim[0] << "; and " << problo[2] << " z < " << m_gran_lim[1] << std::endl;
#elif (AMREX_SPACEDIM == 2)
    const Real H = probhi[1] - problo[1];
    amrex::Print() << "granular region: " << problo[0] << " < x < " << m_gran_lim[0] << "; and " << problo[1] << " y < " << m_gran_lim[1] << std::endl;
#endif
    
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
        Real y = problo[1] + (j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + (k+0.5)*dx[2];

        if ((x < m_gran_lim[0]) and (z < m_gran_lim[1])) {
            density(i,j,k) = rho_2;
            tracer(i,j,k) = 1.0;
        }
        else {
            density(i,j,k) = rho_1;
            tracer(i,j,k) = 0.0;
        }
        gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
        gp(i,j,k,1) = 0.0; 
        gp(i,j,k,2) = m_gravity[2] * density(i,j,k);

        gp0(i,j,k,0) = 0.0;
        gp0(i,j,k,1) = 0.0;
        gp0(i,j,k,2) = 0.0;

        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
        if ((x < m_gran_lim[0]) and (y < m_gran_lim[1])) {
            density(i,j,k) = rho_2;
            tracer(i,j,k) = 1.0;
        }
        else {
            density(i,j,k) = rho_1;
            tracer(i,j,k) = 0.0;
        }
        gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
        gp(i,j,k,1) = m_gravity[1] * density(i,j,k);
        
        gp0(i,j,k,0) = 0.0;
        gp0(i,j,k,1) = 0.0;

        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
#endif

        density0(i,j,k) = density(i,j,k);
    });

    amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = problo[0] + i*dx[0];
        Real y = problo[1] + j*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + k*dx[2];
#endif
#if (AMREX_SPACEDIM == 3)
        Real p_at_split = m_p_top_surface + (H-m_gran_lim[1])*rho_1*std::abs(m_gravity[2]); // pressure at interface
        if (x <  m_gran_lim[0]) { // contains granular interface
            if (z > m_gran_lim[1]) {
                p0(i,j,k) = m_p_top_surface + (H-z)*rho_1*std::abs(m_gravity[2]);
            }
            else {
                p0(i,j,k) = p_at_split + (m_gran_lim[1]-z)*rho_2*std::abs(m_gravity[2]);
            }
        }
        else { // all lighter fluid
            p0(i,j,k) = m_p_top_surface + (H-z)*rho_1*std::abs(m_gravity[2]);
        }
#elif (AMREX_SPACEDIM == 2)
        Real p_at_split = m_p_top_surface + (H-m_gran_lim[1])*rho_1*std::abs(m_gravity[1]); // pressure at interface
        if (x <  m_gran_lim[0]) { // contains granular interface
            if (y > m_gran_lim[1]) {
                p0(i,j,k) = m_p_top_surface + (H-y)*rho_1*std::abs(m_gravity[1]);
            }
            else {
                p0(i,j,k) = p_at_split + (m_gran_lim[1]-y)*rho_2*std::abs(m_gravity[1]);
            }
        }
        else { // all lighter fluid
            p0(i,j,k) = m_p_top_surface + (H-y)*rho_1*std::abs(m_gravity[1]);
        }
#endif
        p_visc(i,j,k) = p0(i,j,k);
        p_nd(i,j,k) = p0(i,j,k);
    });
}


void incflo::inclined_channel (Box const& vbx, Box const& /*gbx*/, Box const& nbx,
                               Array4<Real> const& vel,
                               Array4<Real> const& density,
                               Array4<Real> const& tracer,
                               Array4<Real> const& gp,
                               Box const& /*domain*/,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               GpuArray<Real, AMREX_SPACEDIM> const& problo,
                               GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real rho_1, rho_2, rho;
    Real mu_1, mu_2, nu_1, nu_2, nu;

#if (AMREX_SPACEDIM == 3)
    const Real split = 0.5*(problo[2] + probhi[2]);
    const Real H = probhi[2] - problo[2];
#elif (AMREX_SPACEDIM == 2)
    const Real split = 0.5*(problo[1] + probhi[1]);
    const Real H = probhi[1] - problo[1];
#endif

    if (m_do_vof) {
        rho_1 = m_fluid_vof[0].rho;
        rho_2 = m_fluid_vof[1].rho;
        mu_1 = m_fluid_vof[0].mu;
        mu_2 = m_fluid_vof[1].mu;
        nu_1 = m_fluid_vof[0].mu/rho_1; // kinematic viscosity
        nu_2 = m_fluid_vof[1].mu/rho_2; // kinematic viscosity
        nu = std::min(nu_1,nu_2); // kinematic viscosity of top fluid
        amrex::Print() << "rho1 and rho2 during incline channel setup: " << rho_1 << " " << rho_2 << std::endl;
    }
    else {
        rho = m_ro_0;
        nu = m_fluid.mu/m_fluid.rho; // kinematic viscosity
        amrex::Print() << "rho during single_fluid incline channel setup: " << rho << std::endl;
    }

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = problo[1] + (j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + (k+0.5)*dx[2];
#endif
        if (m_do_vof) {
            if (m_smoothing_width < 0.0) { // discontinuous transition
#if (AMREX_SPACEDIM == 3)
                if (z > split) {
                    density(i,j,k) = rho_1;
                    tracer(i,j,k) = 0.0;
                }
                else {
                    density(i,j,k) = rho_2;
                    tracer(i,j,k) = 1.0;
                }
#endif
#if (AMREX_SPACEDIM == 2)
                if (y > split) {
                    density(i,j,k) = rho_1;
                    tracer(i,j,k) = 0.0;
                }
                else {
                    density(i,j,k) = rho_2;
                    tracer(i,j,k) = 1.0;
                }
#endif
            }
            else { // smoothed interface
#if (AMREX_SPACEDIM == 3)
                Real z_rel = problo[2] + (k+0.5)*dx[2] - split;
                Real smoother = 0.5*std::tanh(z_rel/(m_smoothing_width*dx[2]))+0.5; //goes from 0 to 1
#elif (AMREX_SPACEDIM == 2)
                Real y_rel = problo[1] + (j+0.5)*dx[1] - split;
                Real smoother = 0.5*std::tanh(y_rel/(m_smoothing_width*dx[1]))+0.5; //goes from 0 to 1
#endif
                tracer(i,j,k) = 1.0 - smoother;
                density(i,j,k) = rho_1*smoother + rho_2*(1.0-smoother);
            
            }

#if (AMREX_SPACEDIM == 3)
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);          
            gp(i,j,k,1) = 0.0;
            gp(i,j,k,2) = m_gravity[2] * density(i,j,k);

            // two-layer newtonian fluid solution
            Real halfH = 0.5*H;
            Real A = ((-m_gravity[0] * halfH)/(2.0*(mu_1+mu_2))) * (rho_1*mu_2 + rho_2*mu_1 + 2.0*rho_2*mu_2);
            Real C = ((-m_gravity[0] * halfH)/(2.0*(mu_1+mu_2))) * (3.0*rho_1*mu_2 - rho_2*mu_1 + 2.0*rho_1*mu_1);
            if (z<split) {
                vel(i,j,k,0) = (rho_2*m_gravity[0]/(2.0*mu_2))*z*z + A*z/mu_2;
            }
            else {
                vel(i,j,k,0) = (rho_1*m_gravity[0]/(2.0*mu_1))*(z*z - H*H) + C*(z-H)/mu_1;
            }
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
            gp(i,j,k,1) = m_gravity[1] * density(i,j,k);
            
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;

            // two-layer newtonian fluid solution
            Real halfH = 0.5*H;
            Real A = ((-m_gravity[0] * halfH)/(2.0*(mu_1+mu_2))) * (rho_1*mu_2 + rho_2*mu_1 + 2.0*rho_2*mu_2);
            Real C = ((-m_gravity[0] * halfH)/(2.0*(mu_1+mu_2))) * (3.0*rho_1*mu_2 - rho_2*mu_1 + 2.0*rho_1*mu_1);
            if (y<split) {
                vel(i,j,k,0) = (rho_2*m_gravity[0]/(2.0*mu_2))*y*y + A*y/mu_2;
            }
            else {
                vel(i,j,k,0) = (rho_1*m_gravity[0]/(2.0*mu_1))*(y*y - H*H) + C*(y-H)/mu_1;
            }
            vel(i,j,k,1) = 0.0;

#endif
        }
        else {
            density(i,j,k) = rho;
            tracer(i,j,k) = 1.0;

#if (AMREX_SPACEDIM == 3)
            gp(i,j,k,0) = 0.0;
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);          
            gp(i,j,k,1) = 0.0;
            gp(i,j,k,2) = m_gravity[2] * density(i,j,k);
            
            vel(i,j,k,0) = 0.0;
//            vel(i,j,k,0) = (m_gravity[0]*H/nu)*z - (0.5*m_gravity[0]/nu)*z*z;
//            vel(i,j,k,0) = (m_gravity[0]*H/nu)*z;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
//            gp(i,j,k,0) = 0.0;
            gp(i,j,k,1) = m_gravity[1] * density(i,j,k);
            
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,0) = y*(rho*m_gravity[0])*0.5*(probhi[1]-problo[1])/m_fluid.mu;
//            vel(i,j,k,0) = (m_gravity[0]*H/nu)*y - (0.5*m_gravity[0]/nu)*y*y;
//            vel(i,j,k,0) = (m_gravity[0]*H/nu)*y;
            vel(i,j,k,1) = 0.0;
#endif
        }
    });
}

void incflo::inclined_channel_granular (Box const& vbx, Box const& /*gbx*/, Box const& nbx,
                                        Array4<Real> const& vel,
                                        Array4<Real> const& density,
                                        Array4<Real> const& tracer,
                                        Array4<Real> const& p_nd,
                                        Array4<Real> const& p0,
                                        Array4<Real> const& p_visc,
                                        Array4<Real> const& gp,
                                        Array4<Real> const& gp0,
                                        Array4<Real> const& density0,
                                        Box const& /*domain*/,
                                        GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                        GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                        GpuArray<Real, AMREX_SPACEDIM> const& probhi,
                                        bool half_cell)
{
    Real rho_1, rho_2, rho;
    Real nu_1, nu_2, nu;

#if (AMREX_SPACEDIM == 3)
    const Real split = 0.5*(problo[2] + probhi[2]);
    const Real H = probhi[2] - problo[2];
#elif (AMREX_SPACEDIM == 2)
    const Real split = 0.5*(problo[1] + probhi[1]);
    const Real H = probhi[1] - problo[1];
#endif

    if (m_do_vof) {
        rho_1 = m_fluid_vof[0].rho;
        rho_2 = m_fluid_vof[1].rho;
        nu_1 = m_fluid_vof[0].mu/rho_1; // kinematic viscosity
        nu_2 = m_fluid_vof[1].mu/rho_2; // kinematic viscosity
        nu = std::min(nu_1,nu_2); // kinematic viscosity of top fluid
        amrex::Print() << "rho1 and rho2 during incline channel setup: " << rho_1 << " " << rho_2 << std::endl;
    }
    else {
        rho = m_fluid.rho;
        nu = m_fluid.mu/m_fluid.rho; // kinematic viscosity
        amrex::Print() << "rho during single_fluid incline channel setup: " << rho << std::endl;
    }

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = problo[1] + (j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + (k+0.5)*dx[2];
#endif
        if (m_do_vof) {
            if (m_smoothing_width < 0.0) { // discontinuous transition
#if (AMREX_SPACEDIM == 3)
                if (!half_cell) {
                    if (z > split) {
                        density(i,j,k) = rho_1;
                        tracer(i,j,k) = 0.0;
                    }
                    else {
                        density(i,j,k) = rho_2;
                        tracer(i,j,k) = 1.0;
                    }
                }
                else {
                    if (z > split) {
                        if (z-split > 0.5*dx[2]) {
                            density(i,j,k) = rho_1;
                            tracer(i,j,k) = 0.0;
                        }
                        else if (z-split < 0.5*dx[2]) {
                            Real fac = (z+0.5*dx[2]-split)/dx[2];
                            density(i,j,k) = fac*rho_1 + (1.0-fac)*rho_2;
                            tracer(i,j,k) = 1.0 - fac;
                        }
                        else {
                            density(i,j,k) = 0.5*(rho_1 + rho_2);
                            tracer(i,j,k) = 0.5;
                        }
                    }
                    else if (z < split) {
                        if (split-z > 0.5*dx[2]) {
                            density(i,j,k) = rho_2;
                            tracer(i,j,k) = 1.0;
                        }
                        else if (split-z < 0.5*dx[2]) {
                            Real fac = (split-z+0.5*dx[2])/dx[2];
                            density(i,j,k) = fac*rho_2 + (1.0-fac)*rho_1;
                            tracer(i,j,k) = fac;
                        }
                        else {
                            density(i,j,k) = 0.5*(rho_1 + rho_2);
                            tracer(i,j,k) = 0.5;
                        }
                    }
                    else if (z == split) {
                        density(i,j,k) = 0.5*(rho_1+rho_2);
                        tracer(i,j,k) = 0.5;
                    }
                }
#endif
#if (AMREX_SPACEDIM == 2)
                if (!half_cell) {
                    if (y > split) {
                        density(i,j,k) = rho_1;
                        tracer(i,j,k) = 0.0;
                    }
                    else {
                        density(i,j,k) = rho_2;
                        tracer(i,j,k) = 1.0;
                    }
                }
                else {
                    if (y > split) {
                        if (y-split > 0.5*dx[1]) {
                            density(i,j,k) = rho_1;
                            tracer(i,j,k) = 0.0;
                        }
                        else if (y-split < 0.5*dx[1]) {
                            Real fac = (y+0.5*dx[1]-split)/dx[1];
                            density(i,j,k) = fac*rho_1 + (1.0-fac)*rho_2;
                            tracer(i,j,k) = 1.0 - fac;
                        }
                        else {
                            density(i,j,k) = 0.5*(rho_1 + rho_2);
                            tracer(i,j,k) = 0.5;
                        }
                    }
                    else if (y < split) {
                        if (split-y > 0.5*dx[1]) {
                            density(i,j,k) = rho_2;
                            tracer(i,j,k) = 1.0;
                        }
                        else if (split-y < 0.5*dx[1]) {
                            Real fac = (split-y+0.5*dx[1])/dx[1];
                            density(i,j,k) = fac*rho_2 + (1.0-fac)*rho_1;
                            tracer(i,j,k) = fac;
                        }
                        else {
                            density(i,j,k) = 0.5*(rho_1 + rho_2);
                            tracer(i,j,k) = 0.5;
                        }
                    }
                    else if (y == split) {
                        density(i,j,k) = 0.5*(rho_1+rho_2);
                        tracer(i,j,k) = 0.5;
                    }
                }
#endif
                density0(i,j,k) = density(i,j,k);
            }
            else { // smoothed interface
#if (AMREX_SPACEDIM == 3)
                Real z_rel = problo[2] + (k+0.5)*dx[2] - split;
                Real smoother = 0.5*std::tanh(z_rel/(m_smoothing_width*dx[2]))+0.5; //goes from 0 to 1
#elif (AMREX_SPACEDIM == 2)
                Real y_rel = problo[1] + (j+0.5)*dx[1] - split;
                Real smoother = 0.5*std::tanh(y_rel/(m_smoothing_width*dx[1]))+0.5; //goes from 0 to 1
#endif
                tracer(i,j,k) = 1.0 - smoother;
                density(i,j,k) = rho_1*smoother + rho_2*(1.0-smoother);
            
                density0(i,j,k) = density(i,j,k);
            }

#if (AMREX_SPACEDIM == 3)
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);          
            gp(i,j,k,1) = 0.0;
            gp(i,j,k,2) = m_gravity[2] * density(i,j,k);

            gp0(i,j,k,0) = 0.0;
            gp0(i,j,k,1) = 0.0;
            gp0(i,j,k,2) = 0.0;
            
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);
            gp(i,j,k,1) = m_gravity[1] * density(i,j,k);
            
            gp0(i,j,k,0) = 0.0;
            gp0(i,j,k,1) = 0.0;
  
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
#endif
        }
        else {
            density(i,j,k) = rho;
            tracer(i,j,k) = 1.0;

#if (AMREX_SPACEDIM == 3)
            gp(i,j,k,0) = 0.0;
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);          
            gp(i,j,k,1) = 0.0;
            gp(i,j,k,2) = m_gravity[2] * density(i,j,k);
            
            gp0(i,j,k,0) = 0.0;
            gp0(i,j,k,1) = 0.0;
            gp0(i,j,k,2) = 0.0;
  
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;
#elif (AMREX_SPACEDIM == 2)
            gp(i,j,k,0) = 0.0;
            gp(i,j,k,0) = m_gravity[0] * density(i,j,k);          
            gp(i,j,k,1) = 0.0;
            gp(i,j,k,1) = m_gravity[1] * density(i,j,k);
            
            gp0(i,j,k,0) = 0.0;
            gp0(i,j,k,1) = 0.0;

            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
#endif
        }
    });

    amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = problo[1] + j*dx[1];
#if (AMREX_SPACEDIM == 3)
        Real z = problo[2] + k*dx[2];
#endif
        if (m_do_vof) {
#if (AMREX_SPACEDIM == 3)
            Real p_at_split = m_p_top_surface + (H-split)*rho_1*std::abs(m_gravity[2]); // pressure at interface
            if (z > split) {
                p0(i,j,k) = m_p_top_surface + (H-z)*rho_1*std::abs(m_gravity[2]);
            }
            else {
                p0(i,j,k) = p_at_split + (split-z)*rho_2*std::abs(m_gravity[2]);
            }
#elif (AMREX_SPACEDIM == 2)
            Real p_at_split = m_p_top_surface + (H-split)*rho_1*std::abs(m_gravity[1]); // pressure at interface
            if (y > split) {
                p0(i,j,k) = m_p_top_surface + (H-y)*rho_1*std::abs(m_gravity[1]);
            }
            else {
                p0(i,j,k) = p_at_split + (split-y)*rho_2*std::abs(m_gravity[1]);
            }
#endif
            p_visc(i,j,k) = p0(i,j,k);
            p_nd(i,j,k) = p0(i,j,k);
        }
        else {
#if (AMREX_SPACEDIM == 3)
            p0(i,j,k) = m_p_top_surface + (H-z)*rho*std::abs(m_gravity[2]);
#elif (AMREX_SPACEDIM == 2)
            p0(i,j,k) = m_p_top_surface + (H-y)*rho*std::abs(m_gravity[1]);
#endif
            p_visc(i,j,k) = p0(i,j,k);
            p_nd(i,j,k) = p0(i,j,k);
        }
    });
}

void incflo::init_tuscan (Box const& vbx, Box const& /*gbx*/,
                          Array4<Real> const& vel,
                          Array4<Real> const& density,
                          Array4<Real> const& tracer,
                          Box const& domain,
                          GpuArray<Real, AMREX_SPACEDIM> const& /*dx*/,
                          GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                          GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/)
{
    int half_num_cells = domain.length(2) / 2;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AMREX_D_TERM(vel(i,j,k,0) = Real(0.0);,
                     vel(i,j,k,1) = Real(0.0);,
                     vel(i,j,k,2) = Real(0.0););
        density(i,j,k) = Real(1.0);
        if (k <= half_num_cells) {
            tracer(i,j,k) = Real(0.0);
        } else {
            tracer(i,j,k) = Real(0.01);
         }
    });
}

void incflo::init_boussinesq_bubble (Box const& vbx, Box const& /*gbx*/,
                                     Array4<Real> const& vel,
                                     Array4<Real> const& density,
                                     Array4<Real> const& tracer,
                                     Box const& /*domain*/,
                                     GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                     GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                     GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/) const
{
    if (111 == m_probtype)
    {
        constexpr Real m_fourth = Real(0.25);
        constexpr Real   m_half = Real(0.50);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
#if (AMREX_SPACEDIM == 3)
            vel(i,j,k,2) = 0.0;
#endif
            density(i,j,k) = 1.0;

            Real x = Real(i+0.5)*dx[0];
            Real y = Real(j+0.5)*dx[1];
#if (AMREX_SPACEDIM == 2)
            Real r = std::sqrt((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5));
#elif  (AMREX_SPACEDIM == 3)
            Real z = Real(k+0.5)*dx[2];
            Real r = std::sqrt((x-m_half)*(x-m_half) + (y-m_fourth)*(y-m_fourth) + (z-m_fourth)*(z-m_fourth));
#endif
            if (r < Real(0.1))
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else if (112 == m_probtype) {
        constexpr Real m_fourth = Real(0.25);
        constexpr Real   m_half = Real(0.50);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vel(i,j,k,0) = Real(0.0);
            vel(i,j,k,1) = Real(0.0);
            vel(i,j,k,2) = Real(0.0);
            density(i,j,k) = Real(1.0);

            Real x = Real(i+0.5)*dx[0];
            Real y = Real(j+0.5)*dx[1];
            Real z = Real(k+0.5)*dx[2];

            Real r = std::sqrt( (x-m_fourth)*(x-m_fourth) + (y-m_half )*(y-m_half)
                               +(z-m_fourth)*(z-m_fourth));

            if(r < Real(0.1))
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    } else if (113 == m_probtype) {
        constexpr Real m_fourth = Real(0.25);
        constexpr Real   m_half = Real(0.50);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;
            density(i,j,k) = 1.0;

            Real x = Real(i+0.5)*dx[0];
            Real y = Real(j+0.5)*dx[1];
            Real z = Real(k+0.5)*dx[2];

            Real r = std::sqrt((x-m_fourth)*(x-m_fourth) + (y-m_fourth)*(y-m_fourth) + (z-m_half)*(z-m_half));

            if(r < Real(0.1))
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    }
#endif
}

void incflo::init_periodic_tracer (Box const& vbx, Box const& /*gbx*/,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& /*density*/,
                                   Array4<Real> const& tracer,
                                   Box const& /*domain*/,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real L = probhi[0]-problo[0];
    Real C = Real(2.0)*Real(3.1415926535897932) / L;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr Real A = Real(1.0);
        Real x = Real(i+0.5)*dx[0];
        Real y = Real(j+0.5)*dx[1];
        Real z = Real(k+0.5)*dx[2];
        AMREX_D_TERM(vel(i,j,k,0) = Real(1.0);,
                     vel(i,j,k,1) = Real(0.1)*(std::sin(C*(x+z) - Real(0.00042)) + Real(1.0)) * std::exp(y);,
                     vel(i,j,k,2) = Real(0.1)*(std::sin(C*(x+y) - Real(0.00042)) + Real(1.0)) * std::exp(z););
        tracer(i,j,k) = A *(std::sin(C*(y+z) - Real(0.00042)) + Real(1.0)) * std::exp(x);
    });
}

void incflo::init_double_shear_layer (Box const& vbx, Box const& /*gbx*/,
                                      Array4<Real> const& vel,
                                      Array4<Real> const& /*density*/,
                                      Array4<Real> const& tracer,
                                      Box const& /*domain*/,
                                      GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                      GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                      GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/) const
{
    static constexpr Real twopi = Real(2.0) * Real(3.1415926535897932);
    if (21 == m_probtype)
    {
        constexpr Real m_fourth = Real(0.25);
        constexpr Real m_half   = Real(0.5);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5) * dx[0];
            Real y = Real(j+0.5) * dx[1];
            vel(i,j,k,0) = std::tanh(Real(30.0)*(Real(0.25)-amrex::Math::abs(y-Real(0.5))));
            vel(i,j,k,1) = Real(0.05)*std::sin(twopi*x);
#if (AMREX_SPACEDIM == 3)
            vel(i,j,k,2) = Real(0.0);
#endif
            Real r = std::sqrt((x-m_half)*(x-m_half) + (y-m_fourth)*(y-m_fourth));
            if (r < Real(0.1))
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else if (22 == m_probtype)
    {
        constexpr Real m_half = Real(0.5);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = Real(j+0.5) * dx[1];
            Real z = Real(k+0.5) * dx[2];
            vel(i,j,k,1) = std::tanh(Real(30.0)*(Real(0.25)-amrex::Math::abs(z-m_half)));
            vel(i,j,k,2) = Real(0.05)*std::sin(twopi*y);
            vel(i,j,k,0) = Real(0.0);

            Real r = std::sqrt((y-m_half)*(y-m_half) + (z-m_half)*(z-m_half));
            if (r < Real(0.1))
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    }
    else if (23 == m_probtype)
    {
        constexpr Real m_half = Real(0.5);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5) * dx[0];
            Real z = Real(k+0.5) * dx[2];
            vel(i,j,k,2) = std::tanh(Real(30.0)*(Real(0.25)-amrex::Math::abs(x-m_half)));
            vel(i,j,k,0) = Real(0.05)*std::sin(twopi*z);
            vel(i,j,k,1) = Real(0.0);

            Real r = std::sqrt((x-m_half)*(x-m_half) + (z-m_half)*(z-m_half));
            if (r < .1)
                tracer(i,j,k,0) = Real(0.0);
            else
                tracer(i,j,k,0) = Real(0.01);
        });
    }
#endif
    else
    {
        amrex::Abort("Unknown double shear layer m_probtype");
    };
}

void incflo::init_plane_poiseuille (Box const& vbx, Box const& /*gbx*/,
                                    Array4<Real> const& vel,
                                    Array4<Real> const& /*density*/,
                                    Array4<Real> const& tracer,
                                    Box const& domain,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /*dx*/,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /*problo*/,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /*probhi*/) const
{
    Real dxinv = Real(1.0) / domain.length(0);
    Real dyinv = Real(1.0) / domain.length(1);
#if (AMREX_SPACEDIM == 3)
    Real dzinv = Real(1.0) / domain.length(2);
#else
    Real dzinv = 0.0;
#endif
    const auto dhi = amrex::ubound(domain);

    if (31 == m_probtype)
    {
        Real u = m_ic_u;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = Real(j+0.5)*dyinv;
            AMREX_D_TERM(vel(i,j,k,0) = Real(6.0) * u * y * (Real(1.0)-y);,
                         vel(i,j,k,1) = Real(0.0);,
                         vel(i,j,k,2) = Real(0.0););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && i <= dhi.x/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 && i <= dhi.x/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 && i <= dhi.x*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (311 == m_probtype)
    {
        Real u = m_ic_u;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            AMREX_D_TERM(vel(i,j,k,0) = Real(6.0) * u * z * (Real(1.0)-z);,
                         vel(i,j,k,1) = Real(0.0);,
                         vel(i,j,k,2) = Real(0.0););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && i <= dhi.x/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 && i <= dhi.x/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 && i <= dhi.x*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (41 == m_probtype)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            AMREX_D_TERM(vel(i,j,k,0) = Real(0.5)*z;,
                         vel(i,j,k,1) = Real(0.0);,
                         vel(i,j,k,2) = Real(0.0););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && i <= dhi.x/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 && i <= dhi.x/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 && i <= dhi.x*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (32 == m_probtype)
    {
        Real v = m_ic_v;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = Real(k+0.5)*dzinv;
            AMREX_D_TERM(vel(i,j,k,0) = Real(0.0);,
                         vel(i,j,k,1) = Real(6.0) * v * z * (Real(1.0)-z);,
                         vel(i,j,k,2) = Real(0.0););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && j <= dhi.y/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 && j <= dhi.y/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 && j <= dhi.y*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (322 == m_probtype)
    {
        Real v = m_ic_v;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = Real(i+0.5)*dxinv;
            AMREX_D_TERM(vel(i,j,k,0) = Real(0.0);,
                         vel(i,j,k,1) = Real(6.0) * v * x * (Real(1.0)-x);,
                         vel(i,j,k,2) = Real(0.0););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && j <= dhi.y/8)   tracer(i,j,k,0) = Real(1.0);
            if (nt > 1 && j <= dhi.y/2)   tracer(i,j,k,1) = Real(2.0);
            if (nt > 2 && j <= dhi.y*3/4) tracer(i,j,k,2) = Real(3.0);
        });
    }
    else if (33 == m_probtype)
    {
#if (AMREX_SPACEDIM == 3)
        Real w = m_ic_w;
#endif
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 3)
            Real x = Real(i+0.5)*dxinv;
#endif
            AMREX_D_TERM(vel(i,j,k,0) = 0.0;,
                         vel(i,j,k,1) = 0.0;,
                         vel(i,j,k,2) = Real(6.0) * w * x * (Real(1.0)-x););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && k <= dhi.z/8)   tracer(i,j,k,0) = Real(1.0);
            if (nt > 1 && k <= dhi.z/2)   tracer(i,j,k,1) = Real(2.0);
            if (nt > 2 && k <= dhi.z*3/4) tracer(i,j,k,2) = Real(3.0);
        });
    }
    else if (333 == m_probtype)
    {
#if (AMREX_SPACEDIM == 3)
        Real w = m_ic_w;
#endif
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 3)
            Real y = Real(j+0.5)*dyinv;
#endif
            AMREX_D_TERM(vel(i,j,k,0) = 0.0;,
                         vel(i,j,k,1) = 0.0;,
                         vel(i,j,k,2) = Real(6.0) * w * y * (Real(1.0)-y););

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 && k <= dhi.z/8)   tracer(i,j,k,0) = Real(1.0);
            if (nt > 1 && k <= dhi.z/2)   tracer(i,j,k,1) = Real(2.0);
            if (nt > 2 && k <= dhi.z*3/4) tracer(i,j,k,2) = Real(3.0);
        });
    }
    else
    {
        amrex::Abort("Unknown plane poiseuille m_probtype");
    };
}
