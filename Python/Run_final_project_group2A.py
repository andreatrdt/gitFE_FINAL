# Final FE Project: Multivariate Levy Pricing
# Group 2A

# Description:
# The project focus on the comparison of 2 different methods of model
# pricing; an exponential Levy framework with multiple parameters and a
# Black framework which takes consideration of the volatility only.
# Moreover, the project requires the pricing of an exotic derivative
# between the S&P500 and the EURO50 markets.

# Marco Maspes
# Andrea Tarditi
# Matteo Torba

# Libraries built-in
import numpy as np
import datetime as dt
import scipy.io as sio
import warnings

from scipy.optimize import minimize, fsolve
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

# Libraries customized
import Dataset as obj_data
import calibration_functions as calib_fun
import forward_functions as fwd_fun
import pricing_functions as price_fun

# Disable all warnings
warnings.filterwarnings('ignore')

# Start run time computation
import time
start_time = time.time()

# Fix random seed for the simulation
np.random.seed(42)

#####################################################################################
####################           START THE CODE CREATION          #####################
#####################################################################################

# load files
date_settlement = dt.datetime(2023, 7, 9)

spxsx5e_returns = sio.loadmat('SPXSX5Ereturns.mat')

# Settings of initial datasets
data = sio.loadmat("OptionData.mat")
data_USA = obj_data.Dataset(data['mkt'])
data_EU = obj_data.Dataset(data['mkt_EU'])

# Removal of the last American object
data_USA._remove_last_element()

# Extract the components from the 'Returns' key
daily_returns = spxsx5e_returns['Returns'][0, 0]['Daily']
weekly_returns = spxsx5e_returns['Returns'][0, 0]['Weekly']
monthly_returns = spxsx5e_returns['Returns'][0, 0]['Monthly']
quarterly_returns = spxsx5e_returns['Returns'][0, 0]['Quarterly']
semiannually_returns = spxsx5e_returns['Returns'][0, 0]['SemiAnnually']
annually_returns = spxsx5e_returns['Returns'][0, 0]['Annually']

# Completing the EU and USA dataset with dates

# EU expiry dates
expiry_dates_EU = [
    '2023-07-21', '2023-08-18', '2023-09-15', '2023-10-20',
    '2023-11-17', '2023-12-15', '2024-01-19', '2024-02-16',
    '2024-03-15', '2024-06-21', '2024-09-20', '2024-12-20',
    '2025-06-20'
]

expiry_dates_EU = [dt.datetime.strptime(date, '%Y-%m-%d') for date in expiry_dates_EU]

data_EU._list_dates = expiry_dates_EU

# USA expiry dates
expiry_dates_USA = [
    '2023-07-21', '2023-08-18', '2023-09-15', '2023-10-20', '2023-11-17', '2023-12-15',
    '2024-01-19', '2024-02-16', '2024-03-15', '2024-04-19', '2024-05-17', '2024-06-21',
    '2024-07-19', '2024-09-20', '2024-12-20', '2025-06-20', '2025-12-19', '2026-12-18',
    '2027-12-17']

expiry_dates_USA = [dt.datetime.strptime(date, '%Y-%m-%d') for date in expiry_dates_USA]

data_USA._list_dates = expiry_dates_USA

# Year frac conventions
conv_ACT360 = 2
conv_ACT365 = 3
conv_30360_EU = 6

# Compute the historical correlation for daily returns
daily_correlation = calib_fun.hist_corr(annually_returns)

#%% #################################################################################
####################           FORWARD COMPUTATION              #####################
#####################################################################################

# Computation of the forward prices
F0_EU, B_EU = fwd_fun.forward_prices(data_EU)
F0_USA, B_USA = fwd_fun. forward_prices(data_USA)

#%% #################################################################################
####################           DATA REFINEMENT                  #####################
#####################################################################################

# Options selection and calibration
data_EU_OTM = calib_fun.OTM_preprocessing(data_EU, F0_EU)
data_USA_OTM = calib_fun.OTM_preprocessing(data_USA, F0_USA)

data_calib_EU = calib_fun.dataset_preprocessing(data_EU_OTM, F0_EU, B_EU, date_settlement)
data_calib_USA = calib_fun.dataset_preprocessing(data_USA_OTM, F0_USA, B_USA, date_settlement)

#%% #################################################################################
####################           CALIBRATION MARGINALS            #####################
#####################################################################################

# Calibration of the model
x0 = np.array([0.3, -0.5, 0.15, 0.3, -0.5, 0.15])

# Lower and upper bounds for the parameters
lb = np.array([0.01, -np.inf, 0, 0.01, -np.inf, 0])
ub = [None] * len(lb)

# Options for the optimization algorithm
options = {'maxiter': 3000, 'disp': True}

# Combine the constraints
nonlinear_constraints = {
    'type': 'eq',
    'fun': lambda x: calib_fun.nonlinconstr(x)[1]
}

# The calibration function should be passed without calling it
def calibration_function(params):
    return calib_fun.calibration(params, data_calib_EU, data_calib_USA, F0_EU, B_EU, F0_USA, B_USA, date_settlement)

# Minimize the calibration function with constraints and bounds
res = minimize(calibration_function, x0, method='SLSQP', bounds=[(lb[i], ub[i]) for i in range(len(lb))],
               constraints=nonlinear_constraints, options=options)

params_marginals = res.x

# Display of the parameters on the console
params_USA = params_marginals[:3]
params_EU = params_marginals[3:]

# Explicit useful params
k1 = params_USA[0]
k2 = params_EU[0]

#%% #################################################################################
################              DISPLAY PARAMETERS                    #################
#####################################################################################

print('---------------------------------')
print('k_USA:')
print(params_USA[0])
print('theta_USA:')
print(params_USA[1])
print('sigma_USA:')
print(params_USA[2])
print('---------------------------------')
print('k_EU:')
print(params_EU[0])
print('theta_EU:')
print(params_EU[1])
print('sigma_EU:')
print(params_EU[2])


#%% #################################################################################
################         CALIBRATION IDIOSYNC & SYSTEMATIC            ###############
#####################################################################################

# Try to solve a system to find nu_z, nu_1 and nu_2
rho_theoretical_matlab = 0.2911

def equations(x):
    nu_1, nu_2, nu_z = x

    eq1 = k1 - (nu_1 * nu_z)/(nu_1 + nu_z)
    eq2 = k2 - (nu_2 * nu_z)/(nu_2 + nu_z)
    eq3 = np.sqrt(k1 * k2)/nu_z - rho_theoretical_matlab

    return [eq1, eq2, eq3]


x0 = 0.2 * np.ones(3)
sol = fsolve(equations, x0)

nu_USA, nu_EU, nu_z = sol

# Compute the solution of the system
sol = calib_fun.marginal_param(params_USA, params_EU, nu_z)

# Explicitate the parameters
a_USA = sol[0]
a_EU = sol[1]
beta_z = sol[2]
gamma_z = sol[3]

# Remaining parameters computation
beta_USA = params_USA[1] - a_USA * beta_z
gamma_USA = np.sqrt(params_USA[2]**2 - a_USA**2 * gamma_z**2)
beta_EU = params_EU[1] - a_EU * beta_z
gamma_EU = np.sqrt(params_EU[2]**2 - a_EU**2 * gamma_z**2)

# Idiosyncratic parameters USA
idiosync_USA = [nu_USA, beta_USA, gamma_USA, a_USA]

# Idiosyncratic parameters EU
idiosync_EU = [nu_EU, beta_EU, gamma_EU, a_EU]

# Systematic parameters
syst_Z = [nu_z , beta_z, gamma_z]


#%% #################################################################################
################              DISPLAY PARAMETERS                    #################
#####################################################################################

print('---------------------------------')
print('nu_USA:')
print(nu_USA)
print('nu_EU:')
print(nu_EU)
print('nu_z:')
print(nu_z)

print('---------------------------------')
print('beta_USA:')
print(beta_USA)
print('beta_EU:')
print(beta_EU)
print('beta_z:')
print(beta_z)

print('---------------------------------')
print('gamma_USA:')
print(gamma_USA)
print('gamma_EU:')
print(gamma_EU)
print('gamma_z:')
print(gamma_z)

print('---------------------------------')
print('a_USA:')
print(a_USA)
print('a_EU:')
print(a_EU)


#%% #################################################################################
################                  CALIBRATION BLACK                 #################
#####################################################################################

# Initialization of the parameters
x0 = 1e-4

# Calibration of sigma EU
sigma_EU = minimize(lambda sigma: calib_fun.blk_calibration(sigma, data_calib_EU, F0_EU, B_EU, date_settlement), x0,
                     bounds=[(0, 1)], method='SLSQP').x

# Calibration of sigma USA
sigma_USA = minimize(lambda sigma: calib_fun.blk_calibration(sigma, data_calib_USA, F0_USA, B_USA, date_settlement), x0,
                     bounds=[(0, 1)], method='SLSQP').x

#%% #################################################################################
################                  PRICING LEVY                      #################
#####################################################################################

# Calculate rates for USA and EU
rate_USA, TTM_USA = price_fun.interp_pricing_params(data_calib_USA._list_dates, B_USA, date_settlement, 1)
rate_EU, TTM_EU = price_fun.interp_pricing_params(data_calib_EU._list_dates, B_EU, date_settlement, 1)

# Extract spot prices
S0_USA = data_calib_USA._spot
S0_EU = data_calib_EU._spot

# Compute discount at 1 year
B0_Levy = np.exp(-rate_USA * TTM_USA)

# Concatenate spot prices
S0_Levy = np.array([S0_USA, S0_EU])

# Concatenate rates
rates_Levy = np.array([rate_USA, rate_EU])

# Perform stock simulation
St_Levy = price_fun.stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0_Levy, \
                                          rates_Levy, TTM_USA)

# Unpack simulation results
St_USA_Levy = St_Levy[:, 0]
St_EU_Levy = St_Levy[:, 1]

# Compute certificate payoff
indicator_Levy = St_EU_Levy < (0.95 * S0_EU)
certificate_payoff_Levy = np.maximum(St_USA_Levy - S0_USA, 0) * indicator_Levy

# Compute mean price and confidence interval
mean_price_Levy = np.mean(B0_Levy * certificate_payoff_Levy)
std_dev_Levy = np.std(B0_Levy * certificate_payoff_Levy, ddof=1)
n = len(certificate_payoff_Levy)
z = 1.96  # 95% confidence interval
margin_of_error = z * (std_dev_Levy / np.sqrt(n))
IC_Levy = (mean_price_Levy - margin_of_error, mean_price_Levy + margin_of_error)


#%% #################################################################################
################                  PRICING BLACK                     #################
#####################################################################################

# Forward prices
F01_USA = S0_USA * np.exp(rate_USA * TTM_USA)
F01_EU = S0_EU * np.exp(rate_EU * TTM_EU)

# Simulation of the underlying stock prices
St_Black = price_fun.stock_simulation_Black(np.array([sigma_USA, sigma_EU]), np.array([F01_USA, F01_EU]),
                                                np.array([rate_USA, rate_EU]), 0.801, TTM_USA)

# Unpack simulation results
St_USA_Black = St_Black[:, 0]
St_EU_Black = St_Black[:, 1]

# Compute certificate payoff
indicator_Black = St_EU_Black < (0.95 * S0_EU)
certificate_payoff_Black = np.maximum(St_USA_Black - S0_USA, 0) * indicator_Black

B0_black = np.exp(-rate_USA * TTM_USA)

# Compute mean price and confidence interval
mean_price_Black = np.mean(B0_black * certificate_payoff_Black)
std_dev_Black = np.std(B0_black * certificate_payoff_Black, ddof=1)
n = len(certificate_payoff_Black)
z = 1.96  # 95% confidence interval
margin_of_error = z * (std_dev_Black / np.sqrt(n))
IC_Black = (mean_price_Black - margin_of_error, mean_price_Black + margin_of_error)

#%% #################################################################################
################                  DISPLAY PRICES                    #################
#####################################################################################

print('---------------------------------')
print('Levy Model')
print('Mean price: ', mean_price_Levy)
print('Confidence interval: ', IC_Levy)
print('---------------------------------')
print('Black Model')
print('Mean price: ', mean_price_Black)
print('Confidence interval: ',IC_Black)

# End run time computation and display it
end_time = time.time()
print('---------------------------------')
print(f"Run time: {end_time - start_time} seconds")