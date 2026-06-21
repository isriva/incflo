target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

if (INCFLO_FFT)
   target_sources(incflo PRIVATE IncfloStructFact.cpp)
endif ()
