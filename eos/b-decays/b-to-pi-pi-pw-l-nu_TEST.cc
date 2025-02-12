/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/observable.hh>
#include <eos/b-decays/b-to-psd-psd-l-nu.hh>
#include <eos/maths/complex.hh>
#include <iostream>

using namespace test;
using namespace eos;
using namespace std;

class BToPiPiPWLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToPiPiPWLeptonNeutrinoTest() :
            TestCase("b_to_pi_pi_pw_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                /*p["CKM::abs(V_ub)"]              =  3.32e-3;
                p["B->pipi::a^F1_2_1_0@HKVT2025"]    = 0.012;
                p["B->pipi::a^F1_2_1_1@HKVT2025"]    = 0.011;
                p["B->pipi::a^f_2_1_0@HKVT2025"]    = 0.012;
                p["B->pipi::a^f_2_1_1@HKVT2025"]    = 0.012;
                p["B->pipi::a^g_2_1_0@HKVT2025"]    = 0.012;
                p["B->pipi::a^g_2_1_1@HKVT2025"]    = 0.012;*/

                /*p["B->pipi::a^F1_1_2_0@HKVT2025"]    =  -0.0091;
                p["B->pipi::a^F1_1_2_1@HKVT2025"]    =  -0.0121;
                p["B->pipi::a^f_1_2_0@HKVT2025"]    = -0.0305;
                p["B->pipi::a^f_1_2_1@HKVT2025"]    = -0.0636;
                p["B->pipi::a^g_1_2_0@HKVT2025"]    = -0.0498;
                p["B->pipi::a^g_1_2_1@HKVT2025"]    = -0.0810;*/
                
                p["B->pipi::a^F1_2_1_0@HKVT2025"] = 0.33;
                p["B->pipi::a^g_2_1_0@HKVT2025"] = 0.33;
                Options oo
                {
                    { "model",        "CKM" },
                    { "form-factors", "HKvT2025" },
                    { "scattering-amplitudes", "HKvT2025" },
                    { "integration-points", "2048" },
                    { "U",            "u"       },
                    { "q",            "u"       },
                    { "l",            "e"       },
                    { "I1",           "1"       },
                    { "I2",           "1"       },
                    { "C",           "+-"       }
                };

                BToPPLeptonNeutrino test(p, oo);

                //const double eps = 1e-5;
                std::cerr << test.saturation_1_p() << std::endl;
                std::cerr << test.saturation_1_m() << std::endl;
                /*std::cerr << test.integrated_branching_ratio(1e-1, 8.0, 0.27914, 0.6) << std::endl;
                std::cerr << test.integrated_branching_ratio(8.0, 25.0, 0.27914, 0.6) << std::endl;
                std::cerr << test.integrated_branching_ratio(1e-1, 4.0, 0.6, 0.9) << std::endl;
                std::cerr << test.integrated_branching_ratio(4.0, 8.0, 0.6, 0.9) << std::endl;
                std::cerr << test.integrated_branching_ratio(8.0, 21.9, 0.6, 0.9) << std::endl;*/
                //p["pipi->pipi::P1_B_0@GMKPRDEY2011"]    = 0.066;
                //std::cerr << test.integrated_branching_ratio(8.0, 21.9, 0.6, 0.9) - 9.18e-5 << std::endl;
                /*std::cerr << test.integrated_branching_ratio(1e-1, 4.0, 0.81, 1.44) - 0.70e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(4.0, 8.0, 0.81, 1.44) - 0.64e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(8.0, 19.2, 0.81, 1.44) - 1.78e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(1e-1, 4.0, 1.44, 2.52) - 0.91e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(4.0, 8.0, 1.44, 2.25) - 0.72e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(8.0, 16.7, 1.44, 2.25) - 0.93e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(1e-1, 4.0, 2.25, 5.279) - 0.79e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(4.0, 14.3, 2.25, 5.279) - 0.56e-5 << std::endl;
                std::cerr << test.integrated_branching_ratio(8.0, 25.0, 0.36, 0.81) << std::endl;
                std::cerr << test.double_differential_branching_ratio(0.1, 0.5) << " " << test.double_differential_branching_ratio(2.0, 0.5) << " "
                          << test.double_differential_branching_ratio(4.0, 0.5) << " " << test.double_differential_branching_ratio(8.0, 0.5)
                          << " " << test.double_differential_branching_ratio(18.0, 0.5) << " " << test.double_differential_branching_ratio(20.0, 0.5) << std::endl;*/

            }
        }
} b_to_pi_pi_pw_l_nu_test;
