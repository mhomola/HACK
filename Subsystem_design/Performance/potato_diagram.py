from Subsystem_design.common_constants import Constants
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics

import numpy as np
import matplotlib.pyplot as plt

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
# Up to line 214, kerosene and hydrogen are added together and shown as fuel. After line 214 they're shown separately (the code is basically the same though, of course)
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

### Define relevant parameters ### test for the A320 (some rough estimations)

# Weight (all values in [kg])
OEW         = 44560  
MTOW        = 73500  
FuelW_ker   = 10000
FuelW_H2    = 0
#
FuelW = FuelW_ker + FuelW_H2
#    
MaxPLW      =   MTOW-OEW-FuelW
Cargo_fd    =   3170
Cargo_af    =   3170

# Wing configuration
xcg0    =   0.25                # [MAC] CG of OEW
MAC     =   4.19                # [m]
b       =   37.57               # [m]
tap     =   0.240               # [-]
sweep   =   27                  # [deg]
x_LEMAC =   17.8                # [m]
x_frontspar = 0.25              # [MAC]
x_rearspar  = 0.2               # [MAC]
SafeMargin  = 0.02              # [MAC]

# Interior configuration
Num_row = 30                    # number of rows [-]
seatpitch = 31 * 0.0254         # [m]
x_firstrow = 6.3                #[m]

# Cargo hold cg
x_cag_fd = 10.5                # [m] cg position of front cargo
x_cag_af = 24.5                # [m] cg position of rear cargo

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

# calculation
y_MAC   = b/2 * (1 + 2*tap)/(3 + 3*tap)
xcg0_m = x_LEMAC + xcg0 * MAC           # xcg in meters
MaxPLW      =   MTOW-OEW-FuelW
xcg_fuel = x_LEMAC + x_frontspar * MAC + 0.5*(x_rearspar -  x_frontspar) * MAC

# CoM of fuel is calculated by assuming the fuel tank is a trapezoid
# It is found that CoM of fuel tank coincides with MAC

### Calculate cg shift ###

## Step 1: adding cargo added
W1_fd = OEW + Cargo_fd
W1_af = OEW + Cargo_af
xcg1_fd = (xcg0_m * OEW + x_cag_fd * Cargo_fd)/W1_fd
xcg1_af = (xcg0_m * OEW + x_cag_af * Cargo_af)/W1_af

W1      = OEW + Cargo_fd + Cargo_af
xcg1    =  (xcg0_m * OEW + x_cag_fd * Cargo_fd + x_cag_af * Cargo_af)/W1

## Step 2: adding window seats
Mass_person = (MaxPLW - Cargo_fd - Cargo_af)/180

# start from the front
W2_win_f = [W1]
xcg_win_f = [xcg1]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W2_win_f.append(W1 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W2_win_f[-2] * xcg_win_f[-1] + 2 * Mass_person * x_seat_j) / (W1 + 2 * j * Mass_person)
    xcg_win_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W2_win_a = np.array(W2_win_f)[::-1]
xcg_win_a = [xcg_win_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row

    xcg = (W2_win_a[i] * xcg_win_a[-1] - 2 * Mass_person * x_seat_j) / (W2_win_a[0] - 2 * j * Mass_person)
    xcg_win_a.append(xcg)                                # cg location with the 2 guys at j-th row

W2 = W1 + Num_row*2*Mass_person
xcg2 = xcg_win_a[0]

## Step 3: adding aisle seats
W3_ais_f = [W2]
xcg_ais_f = [xcg2]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W3_ais_f.append(W2 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W3_ais_f[-2] * xcg_ais_f[-1] + 2 * Mass_person * x_seat_j) / (W2 + 2 * j * Mass_person)
    xcg_ais_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W3_ais_a = np.array(W3_ais_f)[::-1]
xcg_ais_a = [xcg_ais_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row

    xcg = (W3_ais_a[i] * xcg_ais_a[-1] - 2 * Mass_person * x_seat_j) / (W3_ais_a[0] - 2 * j * Mass_person)
    xcg_ais_a.append(xcg)                                # cg location with the 2 guys at j-th row

W3 = W2 + Num_row*2*Mass_person
xcg3 = xcg_ais_a[0]

## Step 4: adding middle seats
W4_mid_f = [W3]
xcg_mid_f = [xcg3]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W4_mid_f.append(W3 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W4_mid_f[-2] * xcg_mid_f[-1] + 2 * Mass_person * x_seat_j) / (W3 + 2 * j * Mass_person)
    xcg_mid_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W4_mid_a = np.array(W4_mid_f)[::-1]
xcg_mid_a = [xcg_mid_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row
    xcg = (W4_mid_a[i] * xcg_mid_a[-1] - 2 * Mass_person * x_seat_j) / (W4_mid_a[i] - 2*Mass_person)
    xcg_mid_a.append(xcg)                                # cg location with the 2 guys at j-th row

W4 = W3 + Num_row*2*Mass_person
xcg4 = xcg_mid_a[0]


## Step 5: adding fuel
W5 = W4 + FuelW
xcg5 = (W4 * xcg4 + FuelW * xcg_fuel)/W5

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

### Plot ###

# Step 1 plot - cargo
x1_list = (np.array([xcg0_m,xcg1_af,xcg1,xcg1_fd,xcg0_m]) - x_LEMAC ) /MAC
W1_list = np.array([OEW,W1_af,W1,W1_fd,OEW])
plt.scatter(x1_list,W1_list,color = 'tab:purple', marker = 'x')
plt.plot(x1_list,W1_list,color = 'tab:purple',linewidth = '1',label = 'Cargo')

# Step 2 plot - pax window seats
x2_list_f = (np.array(xcg_win_f) - x_LEMAC ) /MAC
W2_list_f = np.array(W2_win_f)
plt.scatter(x2_list_f,W2_list_f,color = 'tab:blue', marker = 'x')
plt.plot(x2_list_f,W2_list_f,color = 'tab:blue',linewidth = '1',label = 'Window seats')

x2_list_a = (np.array(xcg_win_a) - x_LEMAC ) /MAC
W2_list_a = np.array(W2_win_a)
plt.scatter(x2_list_a,W2_list_a,color = 'tab:blue', marker = 'x')
plt.plot(x2_list_a,W2_list_a,color = 'tab:blue',linewidth = '1')

#Step 3 plot - pax aisle seats
x3_list_f = (np.array(xcg_ais_f) - x_LEMAC ) /MAC
W3_list_f = np.array(W3_ais_f)
plt.scatter(x3_list_f,W3_list_f,color = 'tab:green', marker = 'x')
plt.plot(x3_list_f,W3_list_f,color = 'tab:green',linewidth = '1',label = 'Aisle seats')

x3_list_a = (np.array(xcg_ais_a) - x_LEMAC ) /MAC
W3_list_a = np.array(W3_ais_a)
plt.scatter(x3_list_a,W3_list_a,color = 'tab:green', marker = 'x')
plt.plot(x3_list_a,W3_list_a,color = 'tab:green',linewidth = '1')

#Step 4 plot - pax middle seats
x4_list_f = (np.array(xcg_mid_f) - x_LEMAC ) /MAC
W4_list_f = np.array(W4_mid_f)
plt.scatter(x4_list_f,W4_list_f,color = 'tab:orange', marker = 'x')
plt.plot(x4_list_f,W4_list_f,color = 'tab:orange',linewidth = '1',label = 'Middle seats')

x4_list_a = (np.array(xcg_mid_a) - x_LEMAC ) /MAC
W4_list_a = np.array(W4_mid_a)
plt.scatter(x4_list_a,W4_list_a,color = 'tab:orange', marker = 'x')
plt.plot(x4_list_a,W4_list_a,color = 'tab:orange',linewidth = '1')

# Step 5 plot - fuel
x5_list = (np.array([xcg4,xcg5]) - x_LEMAC ) /MAC
W5_list = np.array([W4,W5])
plt.scatter(x5_list,W5_list,color = 'tab:red', marker = 'x')
plt.plot(x5_list,W5_list,color = 'tab:red',linewidth = '1',label = 'Fuel')

# Determine extrema of cg
x_min = np.min([min(x1_list), min(x2_list_f),min(x3_list_f),min(x4_list_f),min(x5_list)])
x_max = np.max([max(x1_list), max(x2_list_a),max(x3_list_a),max(x4_list_a),max(x5_list)])

label_sm = str(str(SafeMargin*100)+ ' % Safety margin' )
plt.axvline(x = x_min, color = 'black', linewidth = 1, linestyle = '--',label = label_sm)
plt.axvline(x = x_max, color = 'black', linewidth = 1, linestyle = '--')
plt.axvline(x = x_min-SafeMargin, color = 'black', linewidth = 1, linestyle = '--')
plt.axvline(x = x_max+SafeMargin, color = 'black', linewidth = 1, linestyle = '--')

#plt.xlim([,]) # eventually set x axis limits
plt.grid(axis='y')
plt.xlabel('$x_{cg}$ [MAC]')
plt.ylabel('Mass [kg]')
plt.title('Loading diagram with OEW cg at '+ str(int(xcg0*100)) +' % MAC')

print('Most fd cg:',x_min*100,'% MAC')
print('Most af cg:',x_max*100,'% MAC')
print('(Excluding the safety margin.)')
plt.legend()
plt.show()

print('Original design OEW cg:',xcg0)

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

##### Basically same code as above, just modified to show separately kerosene and hydrogen in the final graph

### Define relevant parameters ### test for the A320 (some rough estimations)

# Weight (all values in [kg])
OEW         = 44560  
MTOW        = 73500  
FuelW_ker   = 8000
FuelW_H2    = 2000
#
FuelW = FuelW_ker + FuelW_H2
#    
MaxPLW      =   MTOW-OEW-FuelW
Cargo_fd    =   3170
Cargo_af    =   3170

# Wing configuration
xcg0    =   0.25                # [MAC] CG of OEW
MAC     =   4.19                # [m]
b       =   37.57               # [m]
tap     =   0.240               # [-]
sweep   =   27                  # [deg]
x_LEMAC =   17.8                # [m]
x_frontspar = 0.25              # [MAC]
x_rearspar  = 0.2               # [MAC]
SafeMargin  = 0.02              # [MAC]

# Interior configuration
Num_row = 30                    # number of rows [-]
seatpitch = 31 * 0.0254         # [m]
x_firstrow = 6.3                #[m]

# Cargo hold cg
x_cag_fd = 10.5                # [m] cg position of front cargo
x_cag_af = 24.5                # [m] cg position of rear cargo

### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

# calculation
y_MAC   = b/2 * (1 + 2*tap)/(3 + 3*tap)
xcg0_m = x_LEMAC + xcg0 * MAC           # xcg in meters
MaxPLW      =   MTOW-OEW-FuelW
xcg_fuel = x_LEMAC + x_frontspar * MAC + 0.5*(x_rearspar -  x_frontspar) * MAC

# CoM of fuel is calculated by assuming the fuel tank is a trapezoid
# It is found that CoM of fuel tank coincides with MAC

### Calculate cg shift ###

## Step 1: adding cargo added
W1_fd = OEW + Cargo_fd
W1_af = OEW + Cargo_af
xcg1_fd = (xcg0_m * OEW + x_cag_fd * Cargo_fd)/W1_fd
xcg1_af = (xcg0_m * OEW + x_cag_af * Cargo_af)/W1_af

W1      = OEW + Cargo_fd + Cargo_af
xcg1    =  (xcg0_m * OEW + x_cag_fd * Cargo_fd + x_cag_af * Cargo_af)/W1

## Step 2: adding window seats
Mass_person = (MaxPLW - Cargo_fd - Cargo_af)/180

# start from the front
W2_win_f = [W1]
xcg_win_f = [xcg1]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W2_win_f.append(W1 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W2_win_f[-2] * xcg_win_f[-1] + 2 * Mass_person * x_seat_j) / (W1 + 2 * j * Mass_person)
    xcg_win_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W2_win_a = np.array(W2_win_f)[::-1]
xcg_win_a = [xcg_win_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row

    xcg = (W2_win_a[i] * xcg_win_a[-1] - 2 * Mass_person * x_seat_j) / (W2_win_a[0] - 2 * j * Mass_person)
    xcg_win_a.append(xcg)                                # cg location with the 2 guys at j-th row

W2 = W1 + Num_row*2*Mass_person
xcg2 = xcg_win_a[0]

## Step 3: adding aisle seats
W3_ais_f = [W2]
xcg_ais_f = [xcg2]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W3_ais_f.append(W2 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W3_ais_f[-2] * xcg_ais_f[-1] + 2 * Mass_person * x_seat_j) / (W2 + 2 * j * Mass_person)
    xcg_ais_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W3_ais_a = np.array(W3_ais_f)[::-1]
xcg_ais_a = [xcg_ais_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row

    xcg = (W3_ais_a[i] * xcg_ais_a[-1] - 2 * Mass_person * x_seat_j) / (W3_ais_a[0] - 2 * j * Mass_person)
    xcg_ais_a.append(xcg)                                # cg location with the 2 guys at j-th row

W3 = W2 + Num_row*2*Mass_person
xcg3 = xcg_ais_a[0]

## Step 4: adding middle seats
W4_mid_f = [W3]
xcg_mid_f = [xcg3]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row, starting from i = 0
    W4_mid_f.append(W3 + 2*j*Mass_person)             # mass with the 2 guys at j-th row
    xcg = (W4_mid_f[-2] * xcg_mid_f[-1] + 2 * Mass_person * x_seat_j) / (W3 + 2 * j * Mass_person)
    xcg_mid_f.append(xcg)                                # cg location with the 2 guys at j-th row

# start from the back
W4_mid_a = np.array(W4_mid_f)[::-1]
xcg_mid_a = [xcg_mid_f[-1]]
for i in range(Num_row):
    j = i + 1                                       # i starts from 0, j starts from 1
    x_seat_j = x_firstrow + i * seatpitch           # cg location of j-th row
    xcg = (W4_mid_a[i] * xcg_mid_a[-1] - 2 * Mass_person * x_seat_j) / (W4_mid_a[i] - 2*Mass_person)
    xcg_mid_a.append(xcg)                                # cg location with the 2 guys at j-th row

W4 = W3 + Num_row*2*Mass_person
xcg4 = xcg_mid_a[0]


## Step 5: adding fuel

# kerosene
W5 = W4 + FuelW_ker
xcg5 = (W4 * xcg4 + FuelW_ker * xcg_fuel)/W5

# hydrogen
W6 = W5 + FuelW_H2
xcg6 = (W5 * xcg5 + FuelW_H2 * xcg_fuel)/W6


### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###

### Plot ###

# Step 1 plot - cargo
x1_list = (np.array([xcg0_m,xcg1_af,xcg1,xcg1_fd,xcg0_m]) - x_LEMAC ) /MAC
W1_list = np.array([OEW,W1_af,W1,W1_fd,OEW])
plt.scatter(x1_list,W1_list,color = 'tab:purple', marker = 'x')
plt.plot(x1_list,W1_list,color = 'tab:purple',linewidth = '1',label = 'Cargo')

# Step 2 plot - pax window seats
x2_list_f = (np.array(xcg_win_f) - x_LEMAC ) /MAC
W2_list_f = np.array(W2_win_f)
plt.scatter(x2_list_f,W2_list_f,color = 'tab:blue', marker = 'x')
plt.plot(x2_list_f,W2_list_f,color = 'tab:blue',linewidth = '1',label = 'Window seats')

x2_list_a = (np.array(xcg_win_a) - x_LEMAC ) /MAC
W2_list_a = np.array(W2_win_a)
plt.scatter(x2_list_a,W2_list_a,color = 'tab:blue', marker = 'x')
plt.plot(x2_list_a,W2_list_a,color = 'tab:blue',linewidth = '1')

#Step 3 plot - pax aisle seats
x3_list_f = (np.array(xcg_ais_f) - x_LEMAC ) /MAC
W3_list_f = np.array(W3_ais_f)
plt.scatter(x3_list_f,W3_list_f,color = 'tab:green', marker = 'x')
plt.plot(x3_list_f,W3_list_f,color = 'tab:green',linewidth = '1',label = 'Aisle seats')

x3_list_a = (np.array(xcg_ais_a) - x_LEMAC ) /MAC
W3_list_a = np.array(W3_ais_a)
plt.scatter(x3_list_a,W3_list_a,color = 'tab:green', marker = 'x')
plt.plot(x3_list_a,W3_list_a,color = 'tab:green',linewidth = '1')

#Step 4 plot - pax middle seats
x4_list_f = (np.array(xcg_mid_f) - x_LEMAC ) /MAC
W4_list_f = np.array(W4_mid_f)
plt.scatter(x4_list_f,W4_list_f,color = 'tab:orange', marker = 'x')
plt.plot(x4_list_f,W4_list_f,color = 'tab:orange',linewidth = '1',label = 'Middle seats')

x4_list_a = (np.array(xcg_mid_a) - x_LEMAC ) /MAC
W4_list_a = np.array(W4_mid_a)
plt.scatter(x4_list_a,W4_list_a,color = 'tab:orange', marker = 'x')
plt.plot(x4_list_a,W4_list_a,color = 'tab:orange',linewidth = '1')

# Step 5 plot - fuel
x5_list = (np.array([xcg4,xcg5]) - x_LEMAC ) /MAC
W5_list = np.array([W4,W5])
plt.scatter(x5_list,W5_list,color = 'tab:red', marker = 'x')
plt.plot(x5_list,W5_list,color = 'tab:red',linewidth = '1',label = 'Kerosene')

x6_list = (np.array([xcg5,xcg6]) - x_LEMAC ) /MAC
W6_list = np.array([W5,W6])
plt.scatter(x6_list,W6_list,color = 'tab:olive', marker = 'x')
plt.plot(x6_list,W6_list,color = 'tab:olive',linewidth = '1',label = 'Hydrogen')

# Determine extrema of cg
x_min = np.min([min(x1_list), min(x2_list_f),min(x3_list_f),min(x4_list_f),min(x5_list),min(x6_list)])
x_max = np.max([max(x1_list), max(x2_list_a),max(x3_list_a),max(x4_list_a),max(x5_list),max(x6_list)])

label_sm = str(str(SafeMargin*100)+ ' % Safety margin' )
plt.axvline(x = x_min, color = 'black', linewidth = 1, linestyle = '--',label = label_sm)
plt.axvline(x = x_max, color = 'black', linewidth = 1, linestyle = '--')
plt.axvline(x = x_min-SafeMargin, color = 'black', linewidth = 1, linestyle = '--')
plt.axvline(x = x_max+SafeMargin, color = 'black', linewidth = 1, linestyle = '--')

#plt.xlim([,]) # eventually set x axis limits
plt.grid(axis='y')
plt.xlabel('$x_{cg}$ [MAC]')
plt.ylabel('Mass [kg]')
plt.title('Loading diagram with OEW cg at '+ str(int(xcg0*100)) +' % MAC')

print('Most fd cg:',x_min*100,'% MAC')
print('Most af cg:',x_max*100,'% MAC')
print('(Excluding the safety margin.)')
plt.legend()
plt.show()

print('Original design OEW cg:',xcg0)