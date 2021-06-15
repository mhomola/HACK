'''Cost analysis from Roskam book chapter 8'''

import numpy as np
from Subsystem_design.common_constants import Constants

#Constants
N_rdte = 5 # number of aircraft produced for the RDTE phase, usually 2-8
F_diff = 2 # difficulty (complexity) of a new aircraft program, F_diff = 1.5 for programs involving moderately agressive
           # used of advanced technology; F_diff = 2.0 for programs involving very aggressive use of advanced technology
F_cad = 0.8 # F_cad = 1.0 for manucatures which are using 'manual' drafting techiques; F_cad = 0.8 for manufactures
            # which are experienced in the use of CAD
R_e_r_1989 =
CEF_thenyear =
CEF_1989 =

C_e_r =     # cost engine
N_e = 2     # number of engines
C_p_r =     # cost propeller
N_p =       # number of propellers
C_avionics_r =    # cost of avionics equipment per aircraft
N_st =      # number of static test aircraft
F_mat = 2.5 # 2-2.5: primary structure is made with 'conventional' composite materials, Li/Al and/or ARAL; 3 for carbon
            # composite airframes
R_t_r =     # tooling labour rate in $/h per manhour
N_r_r = 0.33    # RDTE production rate in units per month. For RDTE a typical rate is 0.33 units per month

F_obs = 1.0 # factor which depends on the importance of having 'low' observables. For commercial aircraft F_obs = 1.0

F_tsf = 0.2 # cost adjustment factor which depends on judgement. F_tsf - 0.0 if no extra facilities are require;
            # F_tsf = 0.2 if extensive test and simulation facilities are required. X-29, ATF and B-2 programs

F_pro_r = 0.1 # for a suggested profit of 10%

F_fin_r = 0.1 # to 0.2, depending on the interest rates which are available

C_RDTE = #C_aed_r + C_dst_r + C_fta_r + C_fto_r + C_tsf_r + C_pro_e + C_fin_r
W_TO =  #[lbs]
V_max = #[keas]


#airframe engineering and design cost: C_aed_r
W_ampr = 10**(0.1936 + 0.8645 * np.log10(W_TO))                                                       #[lbs]
R_e_r_thenyear = R_e_r_1989 * (CEF_thenyear/CEF_1989)
MHR_aed_r = 0.0396 * W_ampr * 0.791 * V_max * 1.526 * N_rdte * 0.183 * F_diff * F_cad
C_aed_r = MHR_aed_r * R_e_r_thenyear

#development support and testing cost: C_dst_r
C_dst_r = 0.008325 * W_ampr * 0.873 * V_max * 1.890 * N_rdte * 0.346 * F_diff
C_e_a_r = (C_e_r * N_e + C_p_r +C_avionics_r)*(N_rdte - N_st)   # cost of engine and avionics as acquired from vendors
                                                                # C_(e+a)_r
MHR_man_r = 28.984 * W_ampr * 0.740 * V_max * 0.543 * N_rdte * 0.524 * F_diff
C_man_r = MHR_man_r * R_m_r
C_mat_r = 37.632*F_mat * W_ampr**(0.689) * V_max**(0.624) * N_rdte**(0.792) * CEF
MHR_tool_r = 4.0127 * W_ampr**(0.764) * V_max**(0.899) * N_rdte**(0.178) * N_r_r**(0.066) * F_diff
C_tool_r = MHR_tool_r * R_t_r
C_qc_r = 0.13 * C_man_r  # quality control cost associated with manufacturing of flight test aircaft
C_fta_r = C_e_a_r + C_man_r + C_mat_r + C_tool_r + C_qc_r

#flight test operations cost: C_fto_r
C_fto_r = 0.001244 * W_ampr**(1.160) * V_max**(1.371) * (N_rdte - N_st)**(1.281) * CEF * F_diff * F_obs

#Test and simulation facilities cost: C_tsf_r
C_tsf_r = F_tsf * C_RDTE

#RDTE profit: C_pro_r
C_pro_r = F_pro_r * C_RDTE

#cost to finance the RDTE phases: C_fin_r
C_fin_r = F_fin_r * C_RDTE

C_ACQ = C_MAN + C_PROC_ACQ = C_MAN + C_PRO
LCC = C_RDTE + C_ACQ + C_OPS + C_DISP

cost_2021 = cost_19xx * (cost_2021/cost_19xx)
