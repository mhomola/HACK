import numpy as np

yrs = 15 #number of years in operation
daily = 6 #daily number of cycles
SF = 1.5 #safety factor
N = yrs * daily * 365 * SF #number of cycles that the tank must withstand
N = 60000
def sigma_fatigue(N):
    """Calculates the allowable load after N cycles. Relation for Alloy 2024-T4"""
    return 855 * N **(-0.109)

ratio = sigma_fatigue(N)/325 #ratio between yield after N cycles and initial yield
print(sigma_fatigue(N))
print(ratio)
