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
     
         FluidVOF_t fluid0;
         FluidVOF_t fluid1;
     
         amrex::ParmParse ppVOF("incflo.vof");

         amrex::Vector<std::string> fluid_model{{"newtonian", "newtonian"}};;
         ppVOF.queryarr("fluid_model", fluid_model);
         if (fluid_model.size() != 2) {
            amrex::Abort("Need 2 incflo.vof.fluid_model");
         }
         amrex::Vector<amrex::Real> rho;
         ppVOF.queryarr("rho", rho);
         if (rho.size() != 2) {
            amrex::Abort("Need 2 incflo.vof.rho");
         }
         else {
            fluid0.rho = rho[0];
            fluid1.rho = rho[1];
         }
         amrex::Vector<amrex::Real> mu;
         mu[0] = 1.0; mu[1] = 1.0;
         ppVOF.queryarr("mu", mu);
         
         amrex::Vector<amrex::Real> diam;
         ppVOF.queryarr("diam", diam);
         amrex::Vector<amrex::Real> n_0;
         ppVOF.queryarr("n_0", n_0);
         amrex::Vector<amrex::Real> tau_0;
         ppVOF.queryarr("tau_0", tau_0);
         amrex::Vector<amrex::Real> papa_reg;
         ppVOF.queryarr("papa_reg", papa_reg);
         amrex::Vector<amrex::Real> eta_0;
         ppVOF.queryarr("eta_0", eta_0);
         amrex::Vector<amrex::Real> alpha_0;
         ppVOF.queryarr("alpha_0", alpha_0);
         amrex::Vector<amrex::Real> A_0;
         ppVOF.queryarr("A_0", A_0);

         amrex::Vector<amrex::Real> mu_1;
         ppVOF.queryarr("mu_1", mu_1);
         amrex::Vector<amrex::Real> n_1;
         ppVOF.queryarr("n_1", n_1);
         amrex::Vector<amrex::Real> tau_1;
         ppVOF.queryarr("tau_1", tau_1);
         amrex::Vector<amrex::Real> papa_reg_1;
         ppVOF.queryarr("papa_reg_1", papa_reg_1);
         amrex::Vector<amrex::Real> eta_1;
         ppVOF.queryarr("eta_1", eta_1);
         amrex::Vector<amrex::Real> alpha_1;
         ppVOF.queryarr("alpha_1", alpha_1);
         amrex::Vector<amrex::Real> A_1;
         ppVOF.queryarr("A_1", A_1);

         amrex::Vector<amrex::Real> mu_2;
         ppVOF.queryarr("mu_2", mu_2);
         amrex::Vector<amrex::Real> n_2;
         ppVOF.queryarr("n_2", n_2);
         amrex::Vector<amrex::Real> tau_2;
         ppVOF.queryarr("tau_2", tau_2);
         amrex::Vector<amrex::Real> papa_reg_2;
         ppVOF.queryarr("papa_reg_2", papa_reg_2);
         amrex::Vector<amrex::Real> eta_2;
         ppVOF.queryarr("eta_2", eta_2);
         amrex::Vector<amrex::Real> alpha_2;
         ppVOF.queryarr("alpha_2", alpha_2);
         amrex::Vector<amrex::Real> A_2;
         ppVOF.queryarr("A_2", A_2);

         // fluid 0
         if(fluid_model[0] == "newtonian")
         {
             fluid0.fluid_model = FluidModel::Newtonian;
             fluid0.mu = mu[0];
             amrex::Print() << "Newtonian fluid0 with"
                            << " mu = " << fluid0.mu << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "powerlaw")
         {
             fluid0.fluid_model = FluidModel::Powerlaw;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mu_1[0] >= 0.0,
                    "viscosity mu_1 must be positive or zero");
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[0] != 1.0,
                     "No point in using power-law rheology with n = 1");
             AMREX_ALWAYS_ASSERT(n_1[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_1[0] != 1.0,
                     "No point in using power-law rheology with n = 1");
             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             fluid0.mu_1 = mu_1[0];
             fluid0.n_1 = n_1[0];
             amrex::Print() << "Power-law fluid0 with"
                            << " mu = " << fluid0.mu
                            << " n_0 = " << fluid0.n_0
                            << " mu_1 = " << fluid0.mu_1
                            << ", n_1 = " << fluid0.n_1 << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "bingham")
         {
             fluid0.fluid_model = FluidModel::Bingham;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using Bingham rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[0] > 0.0,
                        "Papanastasiou regularisation parameter must be positive");
             
             fluid0.mu = mu[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.papa_reg = papa_reg[0];
             amrex::Print() << "Bingham fluid0 with"
                            << " mu = " << fluid0.mu
                            << ", tau_0 = " << fluid0.tau_0
                            << ", papa_reg = " << fluid0.papa_reg << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "hb")
         {
             fluid0.fluid_model = FluidModel::HerschelBulkley;
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mu_1[0] >= 0.0,
                 "viscosities mu_1 must be positive or zero");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[0] != 1.0,
                     "No point in using Herschel-Bulkley rheology with n = 1");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");
             AMREX_ALWAYS_ASSERT(n_1[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_1[0] != 1.0,
                     "No point in using Herschel-Bulkley rheology with n_1 = 1");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_1[0] > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_1 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg_1[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.papa_reg = papa_reg[0];
             fluid0.mu_1 = mu_1[0];
             fluid0.n_1 = n_1[0];
             fluid0.tau_1 = tau_1[0];
             fluid0.papa_reg_1 = papa_reg_1[0];
             amrex::Print() << "Herschel-Bulkley fluid0 with"
                            << " mu = " << fluid0.mu
                            << ", n = " << fluid0.n_0
                            << ", tau_0 = " << fluid0.tau_0
                            << "  papa_reg = " << fluid0.papa_reg
                            << ", mu_1 = " << fluid0.mu_1
                            << ", n_1 = " << fluid0.n_1
                            << ", tau_1 = " << fluid0.tau_1
                            << ", papa_reg_1 = " << fluid0.papa_reg_1 << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "smd")
         {
             fluid0.fluid_model = FluidModel::deSouzaMendesDutra;
             AMREX_ALWAYS_ASSERT(n_0[0] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT(eta_0[0] > 0.0);

             fluid0.mu = mu[0];
             fluid0.n_0 = n_0[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.eta_0 = eta_0[0];
             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << fluid0.mu
                            << ", n = " << fluid0.n_0
                            << ", tau_0 = " << fluid0.tau_0
                            << ", eta_0 = " << fluid0.eta_0 << " and rho = " << fluid0.rho << std::endl;
         }
         else if(fluid_model[0] == "granular")
         {
             fluid0.fluid_model = FluidModel::Granular;
            
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(diam[0] > 0.0,
                        "Particle diameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(A_0[0] > 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha_0[0] >= 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_1[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg_1[0] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(A_1[0] > 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha_1[0] >= 0.0,
                         "Fitting parameter must be positive");

             fluid0.diam = diam[0];
             fluid0.A_0 = A_0[0];
             fluid0.alpha_0 = alpha_0[0];
             fluid0.tau_0 = tau_0[0];
             fluid0.papa_reg = papa_reg[0];
             fluid0.A_1 = A_1[0];
             fluid0.alpha_1 = alpha_1[0];
             fluid0.tau_1 = tau_1[0];
             fluid0.papa_reg_1 = papa_reg_1[0];
            
             amrex::Print() << "Second-order Granular rheology with"
                            << " diameter diam = " << fluid0.diam
                            << ", A_0 = " << fluid0.A_0
                            << " alpha_0 = " << fluid0.alpha_0
                            << " tau_0 = " << fluid0.tau_0
                            << ", papa_reg = " << fluid0.papa_reg
                            << ", A_1 = " << fluid0.A_1
                            << " alpha_1 = " << fluid0.alpha_1
                            << " tau_1 = " << fluid0.tau_1
                            << ", papa_reg_1 = " << fluid0.papa_reg_1 << " and rho = " << fluid0.rho << std::endl;
         }
         else
         {
             amrex::Abort("Unknown fluid_model for fluid0! Choose either newtonian, powerlaw, bingham, hb, smd");
         }
         m_fluid_vof.push_back(fluid0);

         //  fluid 1
         if(fluid_model[1] == "newtonian")
         {
             fluid1.fluid_model = FluidModel::Newtonian;
             fluid1.mu = mu[1];
             amrex::Print() << "Newtonian fluid1 with"
                            << " mu = " << fluid1.mu << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "powerlaw")
         {
             fluid1.fluid_model = FluidModel::Powerlaw;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[1] != 1.0,
                     "No point in using power-law rheology with n = 1");
             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             amrex::Print() << "Power-law fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0 << " and rho = " << fluid1.rho <<  std::endl;
         }
         else if(fluid_model[1] == "bingham")
         {
             fluid1.fluid_model = FluidModel::Bingham;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using Bingham rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[1] > 0.0,
                        "Papanastasiou regularisation parameter must be positive");
             
             fluid1.mu = mu[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.papa_reg = papa_reg[1];
             amrex::Print() << "Bingham fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", tau_0 = " << fluid1.tau_0
                            << ", papa_reg = " << fluid1.papa_reg << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "hb")
         {
             fluid1.fluid_model = FluidModel::HerschelBulkley;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_0[1] != 1.0,
                     "No point in using Herschel-Bulkley rheology with n = 1");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using Herschel-Bulkley rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.papa_reg = papa_reg[1];
             amrex::Print() << "Herschel-Bulkley fluid1 with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0
                            << ", tau_0 = " << fluid1.tau_0
                            << ", papa_reg = " << fluid1.papa_reg << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "smd")
         {
             fluid1.fluid_model = FluidModel::deSouzaMendesDutra;
             AMREX_ALWAYS_ASSERT(n_0[1] > 0.0);
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");
             AMREX_ALWAYS_ASSERT(eta_0[1] > 0.0);

             fluid1.mu = mu[1];
             fluid1.n_0 = n_0[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.eta_0 = eta_0[1];
             amrex::Print() << "de Souza Mendes-Dutra fluid with"
                            << " mu = " << fluid1.mu
                            << ", n = " << fluid1.n_0
                            << ", tau_0 = " << fluid1.tau_0
                            << ", eta_0 = " << fluid1.eta_0 << " and rho = " << fluid1.rho << std::endl;
         }
         else if(fluid_model[1] == "granular")
         {
             fluid1.fluid_model = FluidModel::Granular;
            
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(diam[1] > 0.0,
                        "Particle diameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_0[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(A_0[1] > 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha_0[1] >= 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tau_1[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(papa_reg_1[1] > 0.0,
                     "Papanastasiou regularisation parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(A_1[1] > 0.0,
                         "Fitting parameter must be positive");

             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha_1[1] >= 0.0,
                         "Fitting parameter must be positive");

             fluid1.diam = diam[1];
             fluid1.A_0 = A_0[1];
             fluid1.alpha_0 = alpha_0[1];
             fluid1.tau_0 = tau_0[1];
             fluid1.papa_reg = papa_reg[1];
             fluid1.A_1 = A_1[1];
             fluid1.alpha_1 = alpha_1[1];
             fluid1.tau_1 = tau_1[1];
             fluid1.papa_reg_1 = papa_reg_1[1];
            
             amrex::Print() << "Second-order Granular rheology with"
                            << " diameter diam = " << fluid1.diam
                            << ", A_0 = " << fluid1.A_0
                            << " alpha_0 = " << fluid1.alpha_0
                            << " tau_0 = " << fluid1.tau_0
                            << ", papa_reg = " << fluid1.papa_reg
                            << ", A_1 = " << fluid1.A_1
                            << " alpha_1 = " << fluid1.alpha_1
                            << " tau_1 = " << fluid1.tau_1
                            << ", papa_reg_1 = " << fluid1.papa_reg_1 << " and rho = " << fluid1.rho << std::endl;
         }
         else
         {
             amrex::Abort("Unknown fluid_model for fluid1! Choose either newtonian, powerlaw, bingham, hb, smd");
         }
         m_fluid_vof.push_back(fluid1);

     }

}
