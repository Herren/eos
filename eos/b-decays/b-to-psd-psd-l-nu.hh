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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: B -> D^* l nu
     */
    class BToPPLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToPPLeptonNeutrino>
    {
        public:
            BToPPLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToPPLeptonNeutrino();

            /*!
             * Saturation of the unitarity bounds by three-body modes
             */
            double saturation_0_m() const;
            double saturation_1_m() const;
            double saturation_1_p() const;

            /*!
             * Observables
             */
            double double_differential_decay_width(const double & q2, const double & s) const;
            double double_differential_decay_width_P(const double & q2, const double & s) const;
            double double_differential_decay_width_D(const double & q2, const double & s) const;
            double double_differential_decay_width_P_L(const double & q2, const double & s) const;
            double double_differential_decay_width_D_L(const double & q2, const double & s) const;
            double double_differential_decay_width_P_T(const double & q2, const double & s) const;
            double double_differential_decay_width_D_T(const double & q2, const double & s) const;
            double double_differential_branching_ratio(const double & q2, const double & s) const;
            double double_differential_branching_ratio_P(const double & q2, const double & s) const;
            double double_differential_branching_ratio_D(const double & q2, const double & s) const;
            double double_differential_branching_ratio_P_L(const double & q2, const double & s) const;
            double double_differential_branching_ratio_D_L(const double & q2, const double & s) const;
            double double_differential_branching_ratio_P_T(const double & q2, const double & s) const;
            double double_differential_branching_ratio_D_T(const double & q2, const double & s) const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max, const double & sqrts_min, const double & sqrts_max) const;
            double fully_integrated_branching_ratio() const;
            double fully_integrated_branching_ratio_P() const;
            double fully_integrated_branching_ratio_D() const;
            double integrated_branching_ratio_q2(const double & q2) const;
            double integrated_branching_ratio_sqrts(const double & sqrts) const;
            double integrated_branching_ratio_q2_P(const double & q2) const;
            double integrated_branching_ratio_sqrts_P(const double & sqrts) const;
            double integrated_branching_ratio_q2_D(const double & q2) const;
            double integrated_branching_ratio_sqrts_D(const double & sqrts) const;
            double integrated_branching_ratio_s_P_L(const double & s) const;
            double integrated_branching_ratio_s_P_T(const double & s) const;
            double integrated_branching_ratio_s_D_L(const double & s) const;
            double integrated_branching_ratio_s_D_T(const double & s) const;

            /*!
             * Pseudo-observables
             */
            double re_g_1_1(const double & q2, const double & s) const;
            double re_g_2_1(const double & q2, const double & s) const;
            double re_g_1_2(const double & q2, const double & s) const;
            double re_g_2_2(const double & q2, const double & s) const;
            double im_g_1_1(const double & q2, const double & s) const;
            double im_g_2_1(const double & q2, const double & s) const;
            double im_g_1_2(const double & q2, const double & s) const;
            double im_g_2_2(const double & q2, const double & s) const;

            double re_f_1_1(const double & q2, const double & s) const;
            double re_f_2_1(const double & q2, const double & s) const;
            double re_f_1_2(const double & q2, const double & s) const;
            double re_f_2_2(const double & q2, const double & s) const;
            double im_f_1_1(const double & q2, const double & s) const;
            double im_f_2_1(const double & q2, const double & s) const;
            double im_f_1_2(const double & q2, const double & s) const;
            double im_f_2_2(const double & q2, const double & s) const;

            /*double re_F1_1_1(const double & q2, const double & s) const;
            double re_F1_2_1(const double & q2, const double & s) const;
            double re_F1_1_2(const double & q2, const double & s) const;
            double re_F1_2_2(const double & q2, const double & s) const;
            double im_F1_1_1(const double & q2, const double & s) const;
            double im_F1_2_1(const double & q2, const double & s) const;
            double im_F1_1_2(const double & q2, const double & s) const;
            double im_F1_2_2(const double & q2, const double & s) const;*/

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_nu;
            static const std::string kinematics_description_phi;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif