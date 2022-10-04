target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   DiffusionScalarOp.cpp
   DiffusionScalarOp.H
   # DiffusionTensorOp.cpp
   # DiffusionTensorOp.H
   DiffusionTensorOp1.cpp
   DiffusionTensorOp1.H
   DiffusionTensorOp2.cpp
   DiffusionTensorOp2.H
   incflo_diffusion.cpp
   )
