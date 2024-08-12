/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Florian Herren
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
#include <eos/form-factors/parametric-hkvdwvt2024-impl.hh>

using namespace test;
using namespace eos;

class BToPiPiHKVDWVT2024FormFactorsTest :
    public TestCase
{
    public:
        BToPiPiHKVDWVT2024FormFactorsTest() :
            TestCase("b_to_pipi_hkvdwvt2024_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->pipi::tp_a@HKVDWVT2024"]  = 29.8;
                p["B->pipi::tp_v@HKVDWVT2024"]  = 29.8;
                p["B->pipi::t0@HKVDWVT2024"]    = 0.0;
                p["B->pipi::sin@HKVDWVT2024"]   = 0.81;
                p["B->pipi::s0@HKVDWVT2024"]    = 0.0729;

                p["mass::pi^+@HKVDWVT2024"]     = 0.13957;
                p["mass::B_u@HKVDWVT2024"]      = 5.279;

                HKVDWVT2024FormFactors<BToPiPi2, PToPP2> ff(p, Options{ });

                Diagnostics  diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.0        ,  eps), // y(s  =  4*0.135^2)
                    std::make_pair( -0.009364   ,  eps), // y(s  =  0.1)
                    std::make_pair(  0.0        ,  eps), // z(q2 =  0) - axial
                    std::make_pair(  0.0        ,  eps), // z(q2 =  0) - vector
                    std::make_pair( -0.101852   ,  eps), // z(q2 = 10) - axial
                    std::make_pair( -0.101852   ,  eps), // z(q2 = 10) - vector

                    std::make_pair(  0.430391   ,  eps), // p_0(z = 0.0, s = 0.1)
                    std::make_pair( -0.069127   ,  eps), // p_1(z = 0.0, s = 0.1)
                    std::make_pair(  0.076481   ,  eps), // p_2(z = 0.0, s = 0.1)
                    std::make_pair( -0.082790   ,  eps), // p_3(z = 0.0, s = 0.1)
                    std::make_pair(  0.088232   ,  eps), // p_4(z = 0.0, s = 0.1)
                    std::make_pair( -0.092995   ,  eps), // p_5(z = 0.0, s = 0.1)

                    std::make_pair(  0.430391   ,  eps), // p_0(z = z(q2 = 10, s = 0.1))
                    std::make_pair( -0.113525   ,  eps), // p_1(z = z(q2 = 10, s = 0.1))
                    std::make_pair(  0.089456   ,  eps), // p_2(z = z(q2 = 10, s = 0.1))
                    std::make_pair( -0.093776   ,  eps), // p_3(z = z(q2 = 10, s = 0.1))
                    std::make_pair(  0.100128   ,  eps), // p_4(z = z(q2 = 10, s = 0.1))
                    std::make_pair( -0.106060   ,  eps), // p_5(z = z(q2 = 10, s = 0.1))

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(s = 0.1), 0)
                    std::make_pair( -0.213216   ,  eps), // p_1(y = y(s = 0.1), 0)
                    std::make_pair(  0.045461   ,  eps), // p_2(y = y(s = 0.1), 0)

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(s = 0.1), 0)
                    std::make_pair( -0.209985   ,  eps), // p_1(y = y(s = 0.1), 0)
                    std::make_pair(  0.038999   ,  eps), // p_2(y = y(s = 0.1), 0)

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(s = 0.1), 0)
                    std::make_pair( -0.206883   ,  eps), // p_1(y = y(s = 0.1), 0)
                    std::make_pair(  0.035121   ,  eps), // p_2(y = y(s = 0.1), 0)                    

                    std::make_pair(  0.017914   ,  eps), // phi_g(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)
                    std::make_pair(  0.018162   ,  eps), // phi_g(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)
                    std::make_pair(  0.790116   ,  eps), // phi_g(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)

                    std::make_pair(  0.000785   ,  eps), // phi_f(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)
                    std::make_pair(  0.000815   ,  eps), // phi_f(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)
                    std::make_pair(  0.041125   ,  eps), // phi_f(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)

                    std::make_pair(  0.000025   ,  eps), // phi_F1(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)
                    std::make_pair(  0.000027   ,  eps), // phi_F1(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)
                    std::make_pair(  0.000797   ,  eps), // phi_F1(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)
                    std::make_pair(  0.001781   ,  eps), // phi_F1(z = z(q2 =  4.0), y = y(s = 0.1), l = 0)

                    std::make_pair(  0.008540   ,  eps), // phi_F2(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)
                    std::make_pair(  0.008658   ,  eps), // phi_F2(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)
                    std::make_pair(  0.217464   ,  eps), // phi_F2(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)
                    std::make_pair(  0.000329   ,  eps), // phi_F2(z = z(q2 =  4.0), y = y(s = 0.1), l = 0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                p["B->pipi::a^F2_2_1_0@HKVDWVT2024"]    = 1.0;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 1, true),        0.775601,            eps);
                p["B->pipi::a^F2_2_1_5@HKVDWVT2024"]    = 1.0;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 1, true),        0.484071,            eps);
                p["B->pipi::a^F2_2_1_10@HKVDWVT2024"]    = -0.2;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 1, true),        0.474560,            eps);
                p["B->pipi::a^F2_2_1_3@HKVDWVT2024"]    = 1.0;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 1, true),        1.250161,            eps);
                p["B->pipi::a^F2_1_2_0@HKVDWVT2024"]    = 1.0;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 2, false),        0.654733,            eps);
                p["B->pipi::a^F2_1_2_5@HKVDWVT2024"]    = 1.0;
                TEST_CHECK_NEARLY_EQUAL(ff.unitarity_integrand_0m(0.5, 2, false),        0.411850,            eps);
            }

        }
} b_to_pipi_hkvdwvt2024_form_factors_test;
