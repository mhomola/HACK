import pandas as pd
import numpy as np
from EnergySplit import LHV_hack, MR_h2, MR_ker, ER_h2, ER_ker

class DataFrame:
    def __init__(self):
        self.common_data = \
            pd.DataFrame(data=np.array([['M0', 0.038,0.038,0.1439,0.2495,0.267,0.2864,0.78,0.78,0.24,0.22,0.17,0.01],
                                        ['h', 0, 0,0,0,450,900,11600,11600,900,450,0,0],
                                        ['Thrust', 14000-0.07*792.557, 120000-0.07*33480.94606,120000-0.07*33480.94606,120000-0.07*33480.94606, 46500-0.07*40265.01961,46500-0.07*40265.01961,
                                        19000-0.07-0.07*40572.21543, 19000-0.07-0.07*40572.21543,34000-0.07*127934.76422,34000-0.07*127934.76422, 14000-0.07*12162.0246,5000],
                                        ['A_eff/fan', 3.5, 3.5,2.9, 1.8,1.7,1.6,0.8, 0.8,1.85,2.,2.6,3.65]]),
                         columns=['parameter', 'taxi_out', 'take_off1','take_off2','take_off3', 'climb1','climb2', 'cruise1','cruise2', 'approach1','approach2', 'taxi_in', 'idle'])
        self.neo =\
            pd.DataFrame(data=np.array([[ 'M0',0.038,0.12, 0.497, 0.78, 0.47, 0.038, 0.01],
                                       ['h',0,0,5800,11600,5800,0,0],
                                       ['Thrust', 14000 - 0.07 * 792.557, 120000 - 0.07 * 33480.946, 46500 - 0.07 * 40265.01961,
                                        19000 - 0.07 - 0.07 * 40572.21543, 34000 - 0.07 * 127934.76422,14000 - 0.07 * 12162.0246, 5000],
                                       ['A_eff/fan', 3.5, 3.5, 1, 0.8, 1.05, 3.5, 3.65],# M2 = 0.5, new Mach # 3.5/4.2
                                        # ['A_eff/fan', 3.5, 3, 0.85, 0.7, 0.85, 3.5, 3.65],  # M2 = 0.4, new Mach
                                        ['eta_inlet', 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99],
                                        ['PR_fan', 1.4, 1.4, 1.403, 1.4206, 1.403, 1.4, 1.4],
                                        ['eta_fan', 0.93, 0.93, 0.9197225, 0.909445, 0.9197225, 0.93, 0.93],
                                        ['BPR', 11.1, 11.1, 11.17213, 11.24426, 11.17213, 11.1, 11.1],
                                        ['eta_LPC', 0.92, 0.92, 0.910095, 0.90019, 0.910095, 0.92, 0.92],
                                        ['eta_HPC', 0.9774, 0.9774, 0.966045, 0.95469, 0.966045, 0.9774, 0.9774],
                                        ['PR_LPC', 2, 2,  2.347095, 2.69419,  2.347095, 2, 2],
                                        ['PR_HPC', 11.92999, 11.92999, 10.8339, 9.73784, 10.8339, 11.92999, 11.92999],
                                        ['eta_mech', 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99],
                                        ['eta_cc', 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995],            # that of LEAP-1B
                                        ['PR_cc', 0.94, 0.94, 0.939765, 0.93953, 0.939765, 0.94, 0.94],
                                        ['T04', 1355.153, 1550, 1475, 1300, 1260, 1355.153, 1370],              #1245
                                        ['eta_LPT', 0.95, 0.94, 0.94025, 0.9405, 0.94025, 0.95, 0.94],
                                        ['eta_HPT', 0.95, 0.94, 0.9364, 0.9328, 0.9364, 0.95, 0.94],
                                        ['PR_LPT', 0.157178, 0.157178, 0.141717, 0.12626, 0.141717, 0.157178, 0.157178],
                                        ['PR_HPT', 0.261226, 0.261226, 0.261526, 0.261826, 0.261526, 0.261226, 0.261226],
                                        # ['eta_nozzle', 0.985, 0.98, 0.9809, 0.9818, 0.9809, 0.985, 0.98],
                                        ['eta_nozzle', 0.991, 0.991, 0.9905, 0.99, 0.9905, 0.991, 0.991],
                                        ['PR_noz_core', 0.992, 0.99, 0.98772, 0.985443, 0.98772, 0.992, 0.99],
                                        ['PR_noz_fan', 0.992, 0.99, 0.988722, 0.987444, 0.988722, 0.992, 0.99],
                                        ['mr_h2', 0, 0, 0, 0, 0, 0, 0],
                                        ['mr_ker', 1, 1, 1, 1, 1, 1, 1],
                                        ['ER_h2', 0, 0, 0, 0, 0, 0, 0],
                                        ['ER_ker', 1, 1, 1, 1, 1, 1, 1],
                                        ['LHV', 43, 43, 43, 43, 43, 43, 43],
                                        ['D_fan', 78, 78, 78, 78, 78, 78, 78]]),
                                columns=['parameter', 'taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in', 'idle'])



        self.hack = \
            pd.DataFrame(data=np.array([['M0',0.038,0.12, 0.497, 0.78, 0.47, 0.038, 0.01],
                                        ['h',0,0,5800,11600,5800,0,0],
                                        ['Thrust', 15000, 120000, 46500, 19000, 34000, 15000, 5000],
                                        ['A_eff/fan', 3.5, 3.5, 1, 0.8, 1.05, 3.5, 3.65],
                                        ['eta_inlet', 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995],
                                        ['PR_fan', 1.35, 1.35, 1.353, 1.3706, 1.353, 1.35, 1.35],
                                        ['eta_fan', 0.950452034, 0.950452034, 0.939948517, 0.929445, 0.939948517, 0.950452034, 0.950452034],
                                        ['BPR', 15.60518878, 15.60518878, 15.70659439, 15.808, 15.70659439, 15.60518878, 15.60518878],
                                        ['eta_LPC', 0.930220065, 0.930220065, 0.920205032, 0.91019, 0.920205032, 0.930220065, 0.930220065],
                                        ['eta_HPC', 0.987637878, 0.987637878, 0.976163939, 0.96469, 0.976163939, 0.987637878, 0.987637878],
                                        ['PR_LPC', 1.7, 1.7, 2.0347095, 2.369419, 1.9347095, 1.7, 1.7],
                                        ['PR_HPC', 18.78063, 18.78063, 16.72263672, 14.81599335, 17.58698544, 18.78063, 18.78063],
                                        ['eta_mech', 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995],
                                        ['eta_cc', 0.9995, 0.9995, 0.9995, 0.9995, 0.9995, 0.9995, 0.9995],
                                        ['PR_cc', 0.955007504, 0.955007504, 0.954768752, 0.95453, 0.954768752, 0.955007504, 0.955007504],
                                        # ['T04', 1325, 1655, 1630, 1475, 1500, 1327, 1590], # prediction: 5K/year * 17 years    1800
                                        ['T04', 1476.957, 1677, 1645, 1475, 1355, 1476.957, 1500],      # H2/ker mix
                                        # ['T04', 1460.716, 1677, 1645, 1475, 1355, 1460.716, 1500],      # H2 or kerosene only
                                        ['eta_LPT', 0.959989367, 0.959989367, 0.960244684, 0.9605, 0.960244684, 0.959989367, 0.959989367],
                                        ['eta_HPT', 0.960154374, 0.960154374, 0.956477187, 0.9528, 0.956477187, 0.960154374, 0.960154374],
                                        # NOT USED SO DIDN'T UPDATE
                                        ['PR_LPT', 0.157178, 0.157178, 0.141717, 0.12626, 0.141717, 0.157178, 0.157178],
                                        ['PR_HPT', 0.261226, 0.261226, 0.261526, 0.261826, 0.261526, 0.261226, 0.261226],
                                        # -------------------------
                                        ['eta_nozzle', 0.9931758, 0.9931758, 0.9940879, 0.995, 0.9940879, 0.9931758, 0.9931758],
                                        ['PR_noz_core', 0.994578073, 0.994578073, 0.99228753, 0.99, 0.99228753, 0.994578073, 0.994578073],
                                        ['PR_noz_fan', 0.997575559, 0.997575559, 0.996338162, 0.995, 0.996338162, 0.997575559, 0.997575559],
                                        ['mr_h2', MR_h2[0], MR_h2[1], MR_h2[2], MR_h2[3], MR_h2[4], MR_h2[5], MR_h2[6]],
                                        ['mr_ker', MR_ker[0], MR_ker[1], MR_ker[2], MR_ker[3], MR_ker[4], MR_ker[5], MR_ker[6]],
                                        ['ER_h2', ER_h2[0], ER_h2[1], ER_h2[2], ER_h2[3], ER_h2[4], ER_h2[5], ER_h2[6]],
                                        ['ER_ker', ER_ker[0], ER_ker[1], ER_ker[2], ER_ker[3], ER_ker[4], ER_ker[5], ER_ker[6]],
                                        ['LHV', LHV_hack[0], LHV_hack[1], LHV_hack[2], LHV_hack[3], LHV_hack[4], LHV_hack[5], LHV_hack[6]],
                                        ['D_fan', 81, 81, 81, 81, 81, 81, 81]], dtype="object"),
                         columns=['parameter', 'taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in', 'idle'])

if __name__ == '__main__':
    df_neo = DataFrame().neo
    df_hack = DataFrame().hack
    col_neo = list(df_neo.columns)
    col_hack = list(df_hack.columns)

    N = open('engine_data_neo.txt', 'w')
    for c in range(len(col_neo)):
        N.write(col_neo[c] + '\t')
    N.write('\n')
    for i in range(len(df_neo.values)):
        for j in range(len(df_neo.values[i])):
            N.write(str(df_neo.values[i][j]) + '\t')
        N.write('\n')
    N.close()

    H = open('engine_data_hack.txt', 'w')
    for c in range(len(col_hack)):
        H.write(col_hack[c] + '\t')
    H.write('\n')
    for i in range(len(df_hack.values)):
        for j in range(len(df_hack.values[i])):
            H.write(str(df_hack.values[i][j]) + '\t')
        H.write('\n')
    H.close()
