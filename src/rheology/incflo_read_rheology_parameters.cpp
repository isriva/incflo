#include <incflo.H>

using namespace amrex;

void incflo::ReadRheologyParameters()
{
     amrex::ParmParse pp0("incflo");
     pp0.query("do_vof", m_do_vof);
     
     // Initialize Rheology Parameters for Single Fluid
     if (!m_do_vof) {
         
         // default
         std::string fluid_model_s = "newtonian";
        
         amrex::ParmParse pp("incflo");
         pp.query("fluid_model", fluid_model_s);
         Real temp;
         if (pp.query("mu", temp)) m_fluid.mu = temp;
         else amrex::Abort("need to specify at least incflo.mu viscosity for rheology");
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.mu >= 0.0, "viscosity (mu) must be positive");
        
         // newtonian
         if(fluid_model_s == "newtonian")
         {
            m_fluid.fluid_model = FluidModel::Newtonian;
            amrex::Print() << "Newtonian fluid with"
                           << " mu = " << m_fluid.mu << std::endl;
         }
     
         // second-order power law
         else if(fluid_model_s == "powerlaw")
         {
            m_fluid.fluid_model = FluidModel::Powerlaw;
            
            pp.query("mu_1", m_fluid.mu_1);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.mu_1 >= 0.0,
                    "viscosity mu_1 must be positive or zero");
            pp.query("n", m_fluid.n_0);
            AMREX_ALWAYS_ASSERT(m_fluid.n_0 > 0.0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.n_0 != 1.0,
                    "specify power-law rheology exponent n != 1.0");
            pp.query("n_1", m_fluid.n_1);
            AMREX_ALWAYS_ASSERT(m_fluid.n_1 > 0.0);

            amrex::Print() << "Power-law fluid with"
                           << " mu = " << m_fluid.mu
                           << " mu_1 = " << m_fluid.mu_1
                           << ", n = " << m_fluid.n_0
                           << ", n_1 = " << m_fluid.n_1 << std::endl;
         }

         // Bingham 
         else if(fluid_model_s == "bingham")
         {
            m_fluid.fluid_model = FluidModel::Bingham;
            pp.query("tau_0", m_fluid.tau_0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.tau_0 > 0.0,
                    "No point in using Bingham rheology with tau_0 = 0");

            pp.query("papa_reg", m_fluid.papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");


            amrex::Print() << "Bingham fluid with"
                            << " mu = " << m_fluid.mu
                            << ", tau_0 = " << m_fluid.tau_0
                            << ", papa_reg = " << m_fluid.papa_reg << std::endl;
         }

         // second-order Herschel-Bulkley
         else if(fluid_model_s == "hb")
         {
            m_fluid.fluid_model = FluidModel::HerschelBulkley;
         
            pp.query("mu_1", m_fluid.mu_1);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.mu_1 >= 0.0,
                 "viscosities mu_1 must be positive or zero");

            pp.query("n", m_fluid.n_0);
            AMREX_ALWAYS_ASSERT(m_fluid.n_0 > 0.0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.n_0 != 1.0,
                    "specify power-law rheology exponent n != 1.0");
            pp.query("n_1", m_fluid.n_1);
            AMREX_ALWAYS_ASSERT(m_fluid.n_1 > 0.0);

            pp.query("tau_0", m_fluid.tau_0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.tau_0 > 0.0,
                 "No point in using Herschel-Bulkley rheology with tau_0 = 0");
            pp.query("tau_1", m_fluid.tau_1);

            pp.query("papa_reg", m_fluid.papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");
            pp.query("papa_reg_1", m_fluid.papa_reg_1);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.papa_reg_1 > 0.0,
                    "Papanastasiou regularisation parameter must be positive");

            amrex::Print() << "Herschel-Bulkley second-order fluid with"
                           << " mu = " << m_fluid.mu
                           << " mu_1 = " << m_fluid.mu_1
                           << ", n = " << m_fluid.n_0
                           << ", n_1 = " << m_fluid.n_1
                           << ", tau_0 = " << m_fluid.tau_0
                           << ", tau_1 = " << m_fluid.tau_1
                           << ", papa_reg = " << m_fluid.papa_reg 
                           << ", papa_reg_1 = " << m_fluid.papa_reg_1 << std::endl;
         }

         // de Souza - Mendes - Dutra
         else if(fluid_model_s == "smd")
         {
             m_fluid.fluid_model = FluidModel::deSouzaMendesDutra;
             
             pp.query("n", m_fluid.n_0);
             AMREX_ALWAYS_ASSERT(m_fluid.n_0 > 0.0);

             pp.query("tau_0", m_fluid.tau_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.tau_0 > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");

             pp.query("eta_0", m_fluid.eta_0);
             AMREX_ALWAYS_ASSERT(m_fluid.eta_0 > 0.0);

             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << m_fluid.mu
                            << ", n = " << m_fluid.n_0
                            << ", tau_0 = " << m_fluid.tau_0
                            << ", eta_0 = " << m_fluid.eta_0 << std::endl;
         }

         // second-order Granular rheology
         else if(fluid_model_s == "granular")
         {
             m_fluid.fluid_model = FluidModel::Granular;
            
             pp.query("diam", m_fluid.diam);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.diam > 0.0,
                        "Particle diameter must be positive");

             pp.query("rho", m_fluid.rho);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.rho > 0.0,
                        "Particle material density must be positive");

             //pp.query("p_bg", m_fluid.p_bg);
             //   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.p_bg > 0.0,
             //           "Background pressure must be positive");

             pp.query("tau_0", m_fluid.tau_0);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.tau_0 > 0.0,
                        "Papanastasiou regularisation parameter must be positive");

             pp.query("papa_reg", m_fluid.papa_reg);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.papa_reg > 0.0,
                        "Papanastasiou regularisation parameter must be positive");

             pp.query("A_0", m_fluid.A_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.A_0 > 0.0,
                         "Fitting parameter must be positive");

             pp.query("alpha_0", m_fluid.alpha_0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.alpha_0 >= 0.0,
                         "Fitting parameter must be positive");

             pp.query("tau_1", m_fluid.tau_1);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.tau_1 > 0.0,
                        "Papanastasiou regularisation parameter must be positive");

             pp.query("papa_reg_1", m_fluid.papa_reg_1);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.papa_reg_1 > 0.0,
                        "Papanastasiou regularisation parameter must be positive");

             pp.query("A_1", m_fluid.A_1);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.A_1 > 0.0,
                         "Fitting parameter must be positive");

             pp.query("alpha_1", m_fluid.alpha_1);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid.alpha_1 >= 0.0,
                         "Fitting parameter must be positive");
            
             amrex::Print() << "Second-order Granular rheology with"
                            << " diameter diam = " << m_fluid.diam
                            << ", mat. density rho = " << m_fluid.rho
                            << ", A_0 = " << m_fluid.A_0
                            << " alpha_0 = " << m_fluid.alpha_0
                            << " tau_0 = " << m_fluid.tau_0
                            << ", papa_reg = " << m_fluid.papa_reg
                            << ", A_1 = " << m_fluid.A_1
                            << " alpha_1 = " << m_fluid.alpha_1
                            << " tau_1 = " << m_fluid.tau_1
                            << ", papa_reg_1 = " << m_fluid.papa_reg_1 << std::endl;
         }

         else
         {
             amrex::Abort("Unknown fluid_model!");
         }

     }

     // Initialize Rheology Parameters for VOF
     if (m_do_vof) 
     {
     
         amrex::ParmParse ppVOF("incflo.vof");

         amrex::Vector<std::string> names;
         ppVOF.queryarr("names", names);
         if (names.size() != 2) amrex::Abort("need two fluid names for VOF");

         for (int i=0; i<2; ++i) {

             std::string name = names[i];
             FLUID_t fluid;

             Real temp; 
             std::string temp_s, fluid_model_s;

             {
                 amrex::ParmParse pp("incflo.vof.fluid_model");
                 if (pp.query(name.c_str(),temp_s)) fluid_model_s = temp_s;
                 else amrex::Abort("need fluid_model for fluid: " + name);
             }
             {
                 amrex::ParmParse pp("incflo.vof.rho");
                 if (pp.query(name.c_str(),temp)) fluid.rho = temp;
                 else amrex::Abort("need density (rho) for fluid: " + name);
             }

             // newtonian
             if(fluid_model_s == "newtonian")
             {
                 fluid.fluid_model = FluidModel::Newtonian;
                 fluid.mu = 1.0;
                 {
                    amrex::ParmParse pp("incflo.vof.mu");
                    pp.query(name.c_str(),fluid.mu);
                 }
                 amrex::Print() << "Newtonian fluid " + name + "  with"
                               << " mu = " << fluid.mu << std::endl;
             }
     
             // Powerlaw
             else if(fluid_model_s == "powerlaw")
             {
                 fluid.fluid_model = FluidModel::Powerlaw;
                 {
                    amrex::ParmParse pp("incflo.vof.mu");
                    pp.query(name.c_str(),fluid.mu);
                 }
                 {
                    amrex::ParmParse pp("incflo.vof.n");
                    pp.query(name.c_str(),fluid.n_0);
                    AMREX_ALWAYS_ASSERT(fluid.n_0 > 0.0);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.n_0 != 1.0,
                            "specify power-law rheology exponent n != 1.0");
                 }

                amrex::Print() << "Power-law fluid with"
                               << " mu = " << fluid.mu
                               << ", n = " << fluid.n_0 <<  std::endl;
             }

             // Bingham 
             else if(fluid_model_s == "bingham")
             {
                fluid.fluid_model = FluidModel::Bingham;
                
                {
                   amrex::ParmParse pp("incflo.vof.mu");
                   pp.query(name.c_str(),fluid.mu);
                }
                {
                   amrex::ParmParse pp("incflo.vof.tau_0");
                   pp.query(name.c_str(),fluid.tau_0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.tau_0 > 0.0,
                           "No point in using Bingham rheology with tau_0 = 0");
                }
                {
                   amrex::ParmParse pp("incflo.vof.papa_reg");
                   pp.query(name.c_str(),fluid.papa_reg);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.papa_reg > 0.0,
                           "Papanastasiou regularisation parameter must be positive");
                }

                amrex::Print() << "Bingham fluid with"
                                << " mu = " << fluid.mu
                                << ", tau_0 = " << fluid.tau_0
                                << ", papa_reg = " << fluid.papa_reg << std::endl;
             }

             // second-order Herschel-Bulkley
             else if(fluid_model_s == "hb")
             {
                fluid.fluid_model = FluidModel::HerschelBulkley;
             
                {
                   amrex::ParmParse pp("incflo.vof.mu");
                   pp.query(name.c_str(),fluid.mu);
                }
                {
                   amrex::ParmParse pp("incflo.vof.tau_0");
                   pp.query(name.c_str(),fluid.tau_0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.tau_0 > 0.0,
                           "No point in using Herschel-Bulkley rheology with tau_0 = 0");
                }
                {
                   amrex::ParmParse pp("incflo.vof.n");
                   pp.query(name.c_str(),fluid.n_0);
                   AMREX_ALWAYS_ASSERT(fluid.n_0 > 0.0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.n_0 != 1.0,
                        "No point in using Herschel-Bulkley rheology with n = 1");
                }
                {
                   amrex::ParmParse pp("incflo.vof.papa_reg");
                   pp.query(name.c_str(),fluid.papa_reg);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.papa_reg > 0.0,
                           "Papanastasiou regularisation parameter must be positive");
                }
                {
                   amrex::ParmParse pp("incflo.vof.mu_1");
                   pp.query(name.c_str(),fluid.mu_1);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.mu_1 >= 0.0,
                        "viscosities mu_1 must be positive or zero");
                }
                {
                   amrex::ParmParse pp("incflo.vof.tau_1");
                   pp.query(name.c_str(),fluid.tau_1);
                }
                {
                   amrex::ParmParse pp("incflo.vof.n_1");
                   pp.query(name.c_str(),fluid.n_1);
                   AMREX_ALWAYS_ASSERT(fluid.n_1 > 0.0);
                }
                {
                   amrex::ParmParse pp("incflo.vof.papa_reg_1");
                   pp.query(name.c_str(),fluid.papa_reg_1);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.papa_reg_1 > 0.0,
                           "Papanastasiou regularisation parameter must be positive");
                }

                amrex::Print() << "Herschel-Bulkley second-order fluid with"
                               << " mu = " << fluid.mu
                               << " mu_1 = " << fluid.mu_1
                               << ", n = " << fluid.n_0
                               << ", n_1 = " << fluid.n_1
                               << ", tau_0 = " << fluid.tau_0
                               << ", tau_1 = " << fluid.tau_1
                               << ", papa_reg = " << fluid.papa_reg 
                               << ", papa_reg_1 = " << fluid.papa_reg_1 << std::endl;
             }

             // de Souza - Mendes - Dutra
             else if(fluid_model_s == "smd")
             {
                 fluid.fluid_model = FluidModel::deSouzaMendesDutra;
                 
                {
                   amrex::ParmParse pp("incflo.vof.mu");
                   pp.query(name.c_str(),fluid.mu);
                }
                {
                   amrex::ParmParse pp("incflo.vof.tau_0");
                   pp.query(name.c_str(),fluid.tau_0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.tau_0 > 0.0,
                           "No point in using SMD rheology with tau_0 = 0");
                }
                {
                   amrex::ParmParse pp("incflo.vof.eta_0");
                   pp.query(name.c_str(),fluid.eta_0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.eta_0 > 0.0,
                           "No point in using SMD rheology with eta_0 = 0");
                }
                {
                   amrex::ParmParse pp("incflo.vof.n");
                   pp.query(name.c_str(),fluid.n_0);
                   AMREX_ALWAYS_ASSERT(fluid.n_0 > 0.0);
                   AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.n_0 != 1.0,
                        "No point in using SMD rheology with n = 1");
                }

                amrex::Print() << "de Souza Mendes-Dutra fluid with"
                               << " mu = " << fluid.mu
                               << ", n = " << fluid.n_0
                               << ", tau_0 = " << fluid.tau_0
                               << ", eta_0 = " << fluid.eta_0 << std::endl;
             }

             // second-order Granular rheology
             else if(fluid_model_s == "granular")
             {
                 fluid.fluid_model = FluidModel::Granular;
                
                 {
                     amrex::ParmParse pp("incflo.vof.diam");
                     pp.query(name.c_str(), fluid.diam);
                     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.diam > 0.0,
                            "Particle diameter must be positive");
                 }
                 //{
                 //    amrex::ParmParse pp("incflo.vof.p_bg");
                 //    pp.query(name.c_str(), fluid.p_bg);
                 //    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.p_bg > 0.0,
                 //           "Particle diameter must be positive");
                 //}
                 {
                    amrex::ParmParse pp("incflo.vof.tau_0");
                    pp.query(name.c_str(),fluid.tau_0);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.tau_0 >= 0.0,
                            "can not use granular rheology with tau_0 < 0");
                 }
                 {
                    amrex::ParmParse pp("incflo.vof.papa_reg");
                    pp.query(name.c_str(),fluid.papa_reg);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.papa_reg > 0.0,
                            "Papanastasiou regularisation parameter must be positive");
                 }
                 {
                     amrex::ParmParse pp("incflo.vof.A_0");
                     pp.query(name.c_str(), fluid.A_0);
                     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.A_0 >= 0.0,
                            "can not use granular rheology with A_0 < 0.0");
                 }
                 {
                     amrex::ParmParse pp("incflo.vof.alpha_0");
                     pp.query(name.c_str(), fluid.alpha_0);
                     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.alpha_0 >= 0.0,
                            "can not use granular rheology with alpha_0 < 0.0");
                 }
                 {
                    amrex::ParmParse pp("incflo.vof.tau_1");
                    pp.query(name.c_str(),fluid.tau_1);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.tau_1 >= 0.0,
                            "can not use granular rheology with tau_1 < 0");
                 }
                 {
                    amrex::ParmParse pp("incflo.vof.papa_reg_1");
                    pp.query(name.c_str(),fluid.papa_reg_1);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.papa_reg_1 > 0.0,
                            "Papanastasiou regularisation parameter must be positive");
                 }
                 {
                     amrex::ParmParse pp("incflo.vof.A_1");
                     pp.query(name.c_str(), fluid.A_1);
                     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.A_1 >= 0.0,
                            "can not use granular rheology with A_1 < 0.0");
                 }
                 {
                     amrex::ParmParse pp("incflo.vof.alpha_1");
                     pp.query(name.c_str(), fluid.alpha_1);
                     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid.alpha_1 >= 0.0,
                            "can not use granular rheology with alpha_1 < 0.0");
                 }
                
                 amrex::Print() << "Second-order Granular rheology with"
                                << " diameter diam = " << fluid.diam
                                << ", mat. density rho = " << fluid.rho
                                << ", A_0 = " << fluid.A_0
                                << " alpha_0 = " << fluid.alpha_0
                                << " tau_0 = " << fluid.tau_0
                                << ", papa_reg = " << fluid.papa_reg
                                << ", A_1 = " << fluid.A_1
                                << " alpha_1 = " << fluid.alpha_1
                                << " tau_1 = " << fluid.tau_1
                                << ", papa_reg_1 = " << fluid.papa_reg_1 << std::endl;
             }

             else
             {
                 amrex::Abort("Unknown fluid_model!");
             }
            
             m_fluid_vof.push_back(fluid);
         }
         
     }
}
