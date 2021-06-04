
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
    years_OPR = np.delete(t.years_full, np.where(t.years_full == 2005)[0][0] )
    OPR = np.array( [ 13, 29, 29.5, 30, 38, 40.5, 42.5, 45, 48.75 ] )
    t.predict('OPR', years_OPR, OPR)

    # BYPASS RATIO
    years_BPR = np.insert(t.years_full, np.where(t.years_full == 2020)[0][0], 2017)
    BPR = np.array([1.5, 4.3, 4.7, 5.1, 6.2, 6.2, 7.5, 9.2, 12.7, 11.6])
    t.predict('OPR', years_BPR, BPR)

    # HIGH PRESSURE COMPRESSOR (HPC) EFFICIENCY


