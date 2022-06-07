#ifndef _CCM_PC_K_H_
#define _CCM_PC_K_H_

#include <AMReX.H>
#include <ccm_pc.H>

namespace amrex {

///
/// Compute the area of a polygon given the vertices
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real 
areaOfPolygon(Real x[], Real y[], int n){
   Real area = 0.0;
   int j = n - 1;
   for (int i = 0; i < n; i++)
   {
      area += (x[j] + x[i]) * (y[j] - y[i]);
      j = i;
   }
   area *= 0.5;
   return std::abs(area);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool 
cellPenetrationCheck(Real x_0, Real y_0, Real X_poly[], Real Y_poly[], int n)
{
   int count = 0;
   Real a;
   for (int j = 0; j < n-1; j++)
   {
      a = (X_poly[j+1]-x_0) * (Y_poly[j]-y_0) - (X_poly[j]-x_0) * (Y_poly[j + 1]-y_0);
      if ( a <= 0 ) count++;
   }

   // Now do j = n-1
   a = (X_poly[0]-x_0) * (Y_poly[n-1]-y_0) - (X_poly[n-1]-x_0) * (Y_poly[0]-y_0);
   if ( a <= 0 ) count++;

   return (count == n);
}

///
/// Compute the force of vertex at (x1,y1) on vertex at (x0,y0)
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void 
computeSpringDashpot(Real x_0,Real y_0,Real u_0,Real v_0,Real x_1,Real y_1,Real u_1,Real v_1, 
                     Real cell_x0, Real ke, Real gamma, Real constant_r0, Real& fx, Real& fy)
{
    Real dist_x = x_0 - x_1;
    Real dist_y = y_0 - y_1;
    
    Real rel_u = u_0 - u_1;
    Real rel_v = v_0 - v_1;
    
    Real r2 = dist_x*dist_x + dist_y*dist_y;

    Real dist_mag = sqrt(r2);

    // Scale the compression/expansion by r_0 to remove the length scale
    Real delta = (dist_mag-cell_x0) / constant_r0;

    RealVect normal(0.0);
    Real dist_mag_inv = 1.0/dist_mag;
    normal[0] = dist_x * dist_mag_inv;
    normal[1] = dist_y * dist_mag_inv;

    // Scale (v dot n) by r_0 to remove the length scale
    Real vel_dot_normal = (rel_u * normal[0] + rel_v * normal[1]) / constant_r0;

    RealVect fn(0.0);
    // calculate the normal contact force // Kn is spring constant
    fn[0] = -ke*delta*normal[0];
    fn[1] = -ke*delta*normal[1];

    RealVect fd(0.0);
    // Damping force // gamma is damping coeffiecient
    fd[0] = -gamma*vel_dot_normal*normal[0];
    fd[1] = -gamma*vel_dot_normal*normal[1];

    // Combine spring and dashpot forces
    fx = fn[0]+fd[0]; 
    fy = fn[1]+fd[1];
}

///
/// Compute the force of line segment between mns and pls on vertex at (x_0,y_0)
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void 
computeLeonardJones(Real x_0,Real y_0,Real x_mns,Real y_mns,Real x_pls,Real y_pls,
                    Real& r_dist, Real fnb, Real k_adhesive, Real& fx, Real& fy, 
                    Real lj_length_scale, Real lj_repulsive_cutoff, Real lj_attractive_cutoff, 
                    Real constant_r0, bool allow_attractive)
{
    Real dist_x_u = x_0 - x_pls;
    Real dist_y_u = y_0 - y_pls;
    Real dist_x_v = x_0 - x_mns;
    Real dist_y_v = y_0 - y_mns;
    Real dist_x_w = x_pls - x_mns;
    Real dist_y_w = y_pls - y_mns;

    Real r_u = dist_x_u * dist_x_u + dist_y_u * dist_y_u;
    Real r_v = dist_x_v * dist_x_v + dist_y_v * dist_y_v;
    Real r_w = dist_x_w * dist_x_w + dist_y_w * dist_y_w;

    Real dist_mag_u = sqrt(r_u);
    Real dist_mag_inv_u = 1.0/dist_mag_u;
    Real dist_mag_v = sqrt(r_v);
    Real dist_mag_inv_v = 1.0/dist_mag_v;
    Real dist_mag_w = sqrt(r_w);
    Real dist_mag_inv_w = 1.0/dist_mag_w;

    RealVect normal_u(0.);
    normal_u[0] = dist_x_u * dist_mag_inv_u;
    normal_u[1] = dist_y_u * dist_mag_inv_u;

    RealVect normal_v(0.);
    normal_v[0] = dist_x_v * dist_mag_inv_v;
    normal_v[1] = dist_y_v * dist_mag_inv_v;

    RealVect normal_w(0.);
    normal_w[0] = dist_x_w * dist_mag_inv_w;
    normal_w[1] = dist_y_w * dist_mag_inv_w;

    Real wu = ( normal_w[0] * normal_u[0] + normal_w[1] * normal_u[1] );
    Real wv = ( normal_w[0] * normal_v[0] + normal_w[1] * normal_v[1] );

    RealVect normal(0.0);
    if ( wu > 0 ){
       r_dist = dist_mag_u;
       normal[0] = normal_u[0];
       normal[1] = normal_u[1];
    }
    else if ( wv < 0 ){
      r_dist = dist_mag_v;
      normal[0] = normal_v[0];
      normal[1] = normal_v[1];
    }
    else{
       r_dist = std::abs(dist_x_w*dist_y_v - dist_y_w * dist_x_v)* dist_mag_inv_w;
       normal[0] = dist_x_v/r_dist - (dist_x_v*dist_x_w + dist_y_v * dist_y_w) / 
                  (dist_x_w*dist_x_w + dist_y_w*dist_y_w)*dist_x_w/r_dist;
       normal[1] = dist_y_v/r_dist - (dist_x_v*dist_x_w + dist_y_v * dist_y_w) /
                  (dist_x_w*dist_x_w + dist_y_w*dist_y_w)*dist_y_w/r_dist;
    }

    // Initialize to zero
    fx = 0.;
    fy = 0.;

    AMREX_ALWAYS_ASSERT(lj_length_scale > lj_repulsive_cutoff);

    Real r_inv = lj_length_scale / r_dist; // This is always > 1 in the repulsive region
                                           //     and always < 1 in the attractive region

    if (r_dist < lj_repulsive_cutoff)
    {
        // Force here is always repulsive because rdist < lj_repulsive_cutoff < lj_length_scale
        Real fnb_mag = (fnb)*(std::pow(r_inv, 12) - std::pow(r_inv, 6))*(constant_r0*constant_r0/(r_dist*r_dist));
        
        // Issue an error if mag_force > 1e10
        if ( fnb_mag > 1e10 ){
              amrex::Print() << "fnb_mag: " << fnb_mag << std::endl;
              Abort("\n Particle to particle repulsion force in same cell is too large!!\n");
        }

        // Each particle updates its force (no need for atomics)
        fx = fnb_mag*normal[0];
        fy = fnb_mag*normal[1];

    } else if ( allow_attractive && (r_dist < lj_attractive_cutoff) ) {

        // Force here is always attractive
        Real tan_length_scale = 100.*0.1122779236;
        Real fnb_mag = (k_adhesive) * std::atan(tan_length_scale * (r_dist - lj_repulsive_cutoff)/constant_r0);

        // Issue an error if mag_force > 1e10
        if ( fnb_mag > 1e10 ){
              amrex::Print() << "fnb_mag: " << fnb_mag << std::endl;
              Abort("\n Particle to particle attractive force in same cell is too large!!\n");
        }

        // Each particle updates its force (no need for atomics)
        fx = -fnb_mag*normal[0];
        fy = -fnb_mag*normal[1];
    }
}

///
/// Compute the force of line segment between mns and pls on vertex at (x_0,y_0)
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void 
computeTheta(Real x_0,Real y_0,Real x_mns,Real y_mns,Real x_pls,Real y_pls, Real& theta,
             Real& normal_sum_x, Real& normal_sum_y)
{
                      Real dist_x_u = x_0 - x_pls;
                      Real dist_y_u = y_0 - y_pls;
                      Real dist_x_v = x_0 - x_mns;
                      Real dist_y_v = y_0 - y_mns;

                      Real dist_mag_u = sqrt(dist_x_u * dist_x_u + dist_y_u * dist_y_u);
                      Real dist_mag_v = sqrt(dist_x_v * dist_x_v + dist_y_v * dist_y_v);

                      Real dist_mag_inv_u = 1.0/dist_mag_u;
                      Real dist_mag_inv_v = 1.0/dist_mag_v;
    
                      RealVect normal_u(0.);
                      normal_u[0] = dist_x_u * dist_mag_inv_u;
                      normal_u[1] = dist_y_u * dist_mag_inv_u;

                      RealVect normal_v(0.);
                      normal_v[0] = dist_x_v * dist_mag_inv_v;
                      normal_v[1] = dist_y_v * dist_mag_inv_v;

                      normal_sum_x = normal_u[0] + normal_v[0];
                      normal_sum_y = normal_u[1] + normal_v[1];

                      // ||u||*||v|| are already included the normal_u and normal_v vectors 
                      Real x = ( normal_u[0] * normal_v[0] + normal_u[1] * normal_v[1] );

                      if (std::abs(x) > 1.0 + 1.e-12)
                         amrex::Abort("Cant take acos of x < -1 or x > 1");
                      else if (x < -1.)
                      {
                         amrex::Print() << "Setting x back to -1; it is now " << x << std::endl;
                         x = -1.;
                      }
                      else if (x > 1.)
                      {
                         amrex::Print() << "Setting x back to 1; it is now " << x << std::endl;
                         x = 1.;
                      }

                      theta = acos(-x);
}
}
#endif
