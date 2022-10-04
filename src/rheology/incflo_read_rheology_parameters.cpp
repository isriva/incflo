#include <incflo.H>

using namespace amrex;

void incflo::ReadRheologyParameters()
{
     amrex::ParmParse pp("incflo");

     std::string fluid_model_s = "newtonian";
     pp.query("fluid_model", fluid_model_s);
     if(fluid_model_s == "newtonian")
     {
         m_fluid_model = FluidModel::Newtonian;
         amrex::Print() << "Newtonian fluid with"
                        << " mu = " << m_mu << std::endl;
     }
     else if(fluid_model_s == "powerlaw")
     {
         m_fluid_model = FluidModel::powerlaw;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                 "No point in using power-law rheology with n = 1");

         amrex::Print() << "Power-law fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0 <<  std::endl;
     }
     else if(fluid_model_s == "powerlaw2")
     {
         m_fluid_model = FluidModel::powerlaw2;
         pp.query("mu_2", m_mu_1);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_mu_2 >= 0.0,
                 "All viscosities (mu and mu_2) must be positive or zero");

         pp.query("n_1", m_n_0);
         pp.query("n_2", m_n_1);

         amrex::Print() << "Second-order power-law fluid with"
                        << " mu = " << m_mu
                        << " mu_2 = " << m_mu_2
                        << ", n1 = " << m_n_0
                        << ", n2 = " << m_n_1 << std::endl;
     }
     else if(fluid_model_s == "bingham")
     {
         m_fluid_model = FluidModel::Bingham;
         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using Bingham rheology with tau_0 = 0");

         pp.query("papa_reg", m_papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");


         amrex::Print() << "Bingham fluid with"
                        << " mu = " << m_mu
                        << ", tau_0 = " << m_tau_0
                        << ", papa_reg = " << m_papa_reg << std::endl;
     }
     else if(fluid_model_s == "hb")
     {
         m_fluid_model = FluidModel::HerschelBulkley;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_n_0 != 1.0,
                 "No point in using Herschel-Bulkley rheology with n = 1");

         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using Herschel-Bulkley rheology with tau_0 = 0");

         pp.query("papa_reg", m_papa_reg);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                 "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Herschel-Bulkley fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0
                        << ", tau_0 = " << m_tau_0
                        << ", papa_reg = " << m_papa_reg << std::endl;
     }
     else if(fluid_model_s == "hb2")
     {
         m_fluid_model = FluidModel::HerschelBulkley2;
         pp.query("mu_2", m_mu_1);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_mu_2 >= 0.0,
                 "All viscosities (mu and mu_2) must be positive or zero");

         pp.query("n_1", m_n_0);
         pp.query("n_2", m_n_1);

         pp.query("tau_1", m_tau_0);
         pp.query("tau_2", m_tau_1);

         pp.query("papa_reg_1", m_papa_reg);
         pp.query("papa_reg_2", m_papa_reg_1);

         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                 "Papanastasiou regularisation parameter must be positive");
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg_1 > 0.0,
                 "Papanastasiou regularisation parameter must be positive");

         amrex::Print() << "Herschel-Bulkley second-order fluid with"
                        << " mu = " << m_mu
                        << " mu_2 = " << m_mu_2
                        << ", n1 = " << m_n_0
                        << ", n2 = " << m_n_1
                        << ", tau_1 = " << m_tau_0
                        << ", tau_2 = " << m_tau_1
                        << ", papa_reg_1 = " << m_papa_reg 
                        << ", papa_reg_2 = " << m_papa_reg_1 << std::endl;
     }
     else if(fluid_model_s == "smd")
     {
         m_fluid_model = FluidModel::deSouzaMendesDutra;
         pp.query("n", m_n_0);
         AMREX_ALWAYS_ASSERT(m_n_0 > 0.0);

         pp.query("tau_0", m_tau_0);
         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_tau_0 > 0.0,
                 "No point in using de Souza Mendes-Dutra rheology with tau_0 = 0");

         pp.query("eta_0", m_eta_0);
         AMREX_ALWAYS_ASSERT(m_eta_0 > 0.0);

         amrex::Print() << "de Souza Mendes-Dutra fluid with"
                        << " mu = " << m_mu
                        << ", n = " << m_n_0
                        << ", tau_0 = " << m_tau_0
                        << ", eta_0 = " << m_eta_0 << std::endl;
     }
     else if(fluid_model_s == "secondorder")
     {
         m_fluid_model = FluidModel::SecondOrder;
         pp.query("mu_2", m_mu_2);
         AMREX_ALWAYS_ASSERT(m_mu_2 > 0.0);

         amrex::Print() << "second order fluid with"
                        << " mu = " << m_mu
                        << ", mu_2 = " << m_mu_2 << std::endl;
     }
     else if(fluid_model_s == "granular")
     {
         m_fluid_model = FluidModel::Granular;
         pp.query("diam", m_diam);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_diam > 0.0,
                    "Particle diameter must be positive");

         pp.query("ro_0", m_ro_0);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ro_0 > 0.0,
                    "Reference density must be positive");

         pp.query("p_bg", m_p_bg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_p_bg > 0.0,
                    "Background pressure must be positive");

         pp.queryarr("delp", m_delp);

         pp.query("papa_reg", m_papa_reg);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_papa_reg > 0.0,
                    "Papanastasiou regularisation parameter must be positive");

         pp.query("mu_1", m_mu_1);
        //     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_mu_1 > 0.0,
        //             "Fitting parameter mu1 must be positive");

         pp.query("A_1", m_A_1);
        //     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_A_1 > 0.0,
        //             "Fitting parameter must be positive");

         pp.query("alpha_1", m_alpha_1);
        //     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_alpha_1 >= 0.0,
        //             "Fitting parameter must be positive");
        
         pp.query("mu_2", m_mu_2);
         pp.query("A_2", m_A_2);
         pp.query("alpha_2", m_alpha_2);

         pp.query("mu_3", m_mu_3);
         pp.query("A_3", m_A_3);
         pp.query("alpha_3", m_alpha_3);

         pp.query("vert_hi", m_vert_hi);
         pp.query("vert_lo", m_vert_lo);
         pp.query("vert_n", m_vert_n);

         amrex::Print() << "Granular stress with"
                        << " mu_1 = " << m_mu_1
                        << ", A_1 = " << m_A_1
                        << ", alpha_1 = " << m_alpha_1
                        << " mu_2 = " << m_mu_2
                        << ", A_2 = " << m_A_2
                        << ", alpha_2 = " << m_alpha_2
                        << " mu_3 = " << m_mu_3
                        << ", A_3 = " << m_A_3
                        << ", alpha_3 = " << m_alpha_3
                        << ", p_bg = " << m_p_bg
                        << ", papa_reg = " << m_papa_reg << std::endl;
     }
     else
     {
         amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd, or granular");
     }
}
