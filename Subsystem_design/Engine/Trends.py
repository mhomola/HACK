import operator

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures


def trend(x,y,degr,string):
    # transforming the data to include another axis
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    x_pred = np.array([2030])
    x_pred = x_pred[:, np.newaxis]

    polynomial_features = PolynomialFeatures(degree=degr)
    x_poly = polynomial_features.fit_transform(x)
    x_poly_pred = polynomial_features.fit_transform(x_pred)

    model = LinearRegression()
    model.fit(x_poly, y)
    y_poly_pred = model.predict(x_poly)
    y_2030 = model.predict(x_poly_pred)
    print('Predicted value of ',string,' in 2030: ',y_2030[0][0])

    rmse = np.sqrt(mean_squared_error(y,y_poly_pred))
    r2 = r2_score(y,y_poly_pred)
    print('RMSE: ',rmse)
    print('R2: ', r2)

    plt.figure()
    plt.xlabel('Year', fontsize=17)
    plt.ylabel(string, fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # sort the values of x before line plot
    sort_axis = operator.itemgetter(0)
    sorted_zip = sorted(zip(x, y_poly_pred), key=sort_axis)
    x, y_poly_pred = zip(*sorted_zip)
    plt.errorbar(x, y_poly_pred, yerr=y_poly_pred - y, fmt='.r', capsize=5)
    plt.plot(x, y_poly_pred, color='k')
    plt.scatter(x, y, s=30)
    plt.show()


    err_arr = np.array([i for i in y_poly_pred - y])
    print('Error array:', err_arr)
    print('Average error:', np.average(np.abs(y_poly_pred - y)))


# DEFINE INPUTS
x_BPR = np.array([1983, 1983, 1984, 1992, 1993, 1996, 2004, 2013, 2014, 2014, 2016, 2016, 2017])
y_BPR = np.array([5.6, 4.46, 5.5, 4.8, 5.6, 5.5, 4.8, 11.1, 12, 12.7, 8.5, 11.1, 11.6])
trend(x_BPR, y_BPR, 3, 'BPR')

x_OPR = np.array([1983, 1983, 1984, 1992, 1993, 1996, 2004, 2013, 2014, 2016, 2016, 2017])
y_OPR = np.array([23, 21.5, 29, 27.5, 31, 22, 25.5, 37.27, 28.5, 42, 33, 38])
trend(x_OPR, y_OPR, 2, 'OPR')

# x_eta_HPC = np.array([1984, 1992, 1993, 1995, 1996, 2013, 2016])
# y_eta_HPC = np.array([0.9125, 0.8205, 0.8725, 0.8, 0.88, 0.95469, 0.92])
# trend(x_eta_HPC, y_eta_HPC, 2, 'eta_HPC')
#
# x_eta_HPT = np.array([1984, 1992, 1993, 1995, 1996, 2013, 2016])
# y_eta_HPT = np.array([0.884, 0.905, 0.9, 0.9125, 0.91, 0.9328, 0.92])
# trend(x_eta_HPT, y_eta_HPT, 2, 'eta_HPT')


x = np.array([6, 15, 24, 33, 41, 52, 59, 66, 73, 81])
y = np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50])

