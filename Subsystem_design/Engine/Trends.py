
""" Improvements in OPR """
import numpy as np
from sklearn.linear_model import LinearRegression

class Trend():
    def __init__(self):
        self.regr = LinearRegression()
        self.years_full = np.arange( 1970, 2035+1, 5 )

    def predict(self, var_string, years, var):
        years = years.reshape(-1,1)
        self.regression = self.regr.fit(years[:len(var)], var)
        self.pred = self.regression.predict(self.years_full.reshape(-1,1))
        self.pred_print = list([['Year', var_string]])
        for i in range(len(self.pred)):
            self.pred_print.append( [ self.years_full[i], self.pred[i] ] )

        print('Prediction for '+var_string+':\n', self.pred_print)

if __name__ == "__main__":
    t = Trend()

    # OVERALL PRESSURE RATIO
    engines_OPR = [ 'JT3D-7', 'CF6-50C', 'PW2040', 'CF6-80C2B1', 'GE90-85B', 'GE90-94B', 'GE90-115B', 'Trent1000-L2', 'Trent XWB-97' ]
    years_OPR = np.delete(t.years_full, np.where(t.years_full == 2005)[0][0] )
    years_OPR = np.delete(years_OPR, np.where(years_OPR == 1990)[0][0])
    OPR = np.array( [ 13, 29, 29.5, 30, 38, 40.5, 42.5, 45, 48.75 ] )
    t.predict('OPR', years_OPR, OPR)

    # BYPASS RATIO
    engines_BPR = [ 'JT3D-7', 'CF6-50C', 'CF6-6D', 'PW2040', 'CF34-3A', 'GE90-85B', 'GE90-94B', 'GEnx-2B67',  'CF6-80C2B1', 'GE90-85B', 'GE90-94B', 'GE90-115B', 'Trent1000-L2', 'Trent XWB-97' ]
    years_BPR = np.insert(t.years_full, np.where(t.years_full == 2020)[0][0], 2017)
    BPR = np.array([1.5, 4.3, 5.9, 5.5, 6.2, 8.3, 8.25, 7.5,    9.2, 12.7, 11.6])
    t.predict('OPR', years_BPR, BPR)

    # HIGH PRESSURE COMPRESSOR (HPC) EFFICIENCY
    years_eta_HPC = np.insert(t.years_full, 0, np.array([1960,1965]))
    years_eta_HPC = np.delete(years_eta_HPC, np.where(years_eta_HPC == 1975)[0][0] )
    years_eta_HPC = np.insert(years_eta_HPC, np.where(years_eta_HPC == 1990)[0][0], 1987)
    years_eta_HPC = np.insert(years_eta_HPC, np.where(years_eta_HPC == 1995)[0][0], 1993)
    eta_HPC = np.array([ 0.754, 0.786, 0.806 ,0.825, 0.821, 0.873, 0.837, 0.919, 0.886, 0.9035, 0.919 ])
    print('years_eta_HPC: ',years_eta_HPC)
    t.predict('Efficiency HPC', years_eta_HPC, eta_HPC)





