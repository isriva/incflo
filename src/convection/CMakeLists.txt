target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_compute_advection_term.cpp
   incflo_compute_MAC_projected_velocities.cpp
   hydro_weno5_edge_state.cpp
   hydro_weno5_extrap_vel_to_faces.cpp
   )
