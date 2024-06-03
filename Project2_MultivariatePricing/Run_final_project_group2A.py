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

# Clearing the workspace
import numpy as np
import scipy.io as sio
import datetime
import warnings
import os
from scipy.optimize import minimize

from functions import dataset_exploration, hist_corr, estimation_discount_factor, \
    synthetic_forward, forward_prices, OTM_preprocessing, yearfrac, blkimpv, \
    blsdelta, dataset_preprocessing, callPriceLewis_pref, RMSE_total,\
    new_calibration, disp_params, disp_marginal_params, interp_pricing_params,\
        nonlinconstr, nonlinconstr_corr,marginal_param, blk_calibration, stock_simulation_Levy,\
            rate_interpolation, stock_simulation_Black, blk_semiclosed, disp_contract_prices



# Disable all warnings
warnings.filterwarnings('ignore')

# Start run time computation
import time
start_time = time.time()

# Fix random seed for the simulation
np.random.seed(42)

# Flags

# flag = 1 for plotting enabled
# flag = 0 for plotting disabled
flag = 0

# save_results = 1 for text saving enabled
# save_results = 0 for text saving disabled
save_results = 1

# Load folders
# Introduced multiple folders for a better division of the scripts
# Folder containing the data files
folder = 'data'

# Append the directory containing the data files to the system path
os.sys.path.append(folder)

# Loading of the matrices
# Loading of the initial matrices and structures
data = sio.loadmat(os.path.join(folder, 'OptionData.mat'))
data_USA = data['mkt']
data_EU = data['mkt_EU']

SP500_EUR50 = sio.loadmat(os.path.join(folder, 'SPXSX5Ereturns.mat'))['Returns']

# # Dates Vector:
# dates_EU = [datetime.datetime.fromordinal(int(date)).date() for date in data_EU['datesExpiry'].flatten()]
# dates_USA = [datetime.datetime.fromordinal(int(date)).date() for date in data_USA['datesExpiry'].flatten()]

# Computation of the historical correlation
rho_historical = hist_corr(SP500_EUR50)

# Year frac conventions
conv_ACT360 = 2
conv_ACT365 = 3
conv_30360_EU = 6

# Plot of initial options
if flag == 1:
    dataset_exploration(data_EU, data_USA, date_settlement)


# Computation of the forward prices
F0_EU, B_EU = forward_prices(data_EU, flag)
F0_USA, B_USA = forward_prices(data_USA, flag)


# Options selection and calibration
data_EU_OTM = OTM_preprocessing(data_EU, F0_EU)
data_USA_OTM = OTM_preprocessing(data_USA, F0_USA)

data_calib_EU = dataset_preprocessing(data_EU_OTM, F0_EU, B_EU, date_settlement, flag)
data_calib_USA = dataset_preprocessing(data_USA_OTM, F0_USA, B_USA, date_settlement, flag)


# Calibration of the model
x0 = np.array([0.3, -0.5, 0.15, 0.3, -0.5, 0.15])
initial_cond = x0

# Linear inequality constraints
A = np.array([[-1, 0, 0, 0, 0, 0],
              [0, 0, -1, 0, 0, 0],
              [0, 0, 0, -1, 0, 0],
              [0, 0, 0, 0, 0, -1]])

b = np.array([0, 0, 0, 0])

# Linear equality constraints
Aeq = None
beq = None

# Lower and upper bounds
lb = np.array([0.01, -np.inf, 0, 0.01, -np.inf, 0])
ub = None

# Options for the visualization
options = {'maxiter': 3000, 'disp': True}

# Calibration
res = minimize(lambda params: new_calibration(params, data_calib_EU, data_calib_USA, F0_EU, B_EU, F0_USA, B_USA, date_settlement),
               x0,
               constraints={'type': 'ineq', 'fun': lambda x: A.dot(x) - b},
               bounds=[(lb[i], None) for i in range(len(lb))],
               method='SLSQP',
               options=options)

params_marginals = res.x

# Display of the parameters on the console
params_USA = params_marginals[:3]
params_EU = params_marginals[3:]

# Explicit useful params
k1 = params_USA[0]
k2 = params_EU[0]

from scipy.optimize import minimize

# Initialization of the parameters
A = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
b = [0, 0, 0]
Aeq = []
beq = []
lb = [0, 0, 0]
ub = []

x0 = [1, 1, 1]

# Calibration of the nu parameters
params = minimize(lambda params: (np.sqrt(params[0] * params[1] / ((params[0] + params[2]) * (params[1] + params[2]))) - rho_historical) ** 2,
                  x0, constraints={'type': 'ineq', 'fun': lambda params: np.array(A) @ params - b},
                  bounds=[(lb[i], ub[i]) for i in range(len(x0))]).x

nu_USA = params[0]
nu_EU = params[1]
nu_z = params[2]


# Display the calibrated parameters
disp_params(params_marginals, initial_cond, save_results)

# Compute rho obtained from the model
rho_model_Levy = np.sqrt(params[0] * params[1] / ((params[0] + params[2]) * (params[1] + params[2])))
error_rho = abs(rho_model_Levy - rho_historical)
print('Error for rho calibration:')
print(error_rho)

# Compute the solution of the system
sol = marginal_param(params_USA, params_EU, nu_z, rho_model_Levy)

# Explicitate the parameters
a_USA = sol.x[0]
a_EU = sol.x[1]
beta_z = sol.x[2]
gamma_z = sol.x[3]

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

# Display of the parameters over the command window
disp_marginal_params(idiosync_USA, idiosync_EU, beta_z, gamma_z, nu_z, save_results)


# BLACK calibration

# Initialization of the parameters
x0 = 1e-4

# Calibration of sigma EU
sigma_EU = minimize(lambda sigma: blk_calibration(sigma, data_calib_EU, F0_EU, B_EU, date_settlement),
                    x0, method='Nelder-Mead').x

# Calibration of sigma USA
sigma_USA = minimize(lambda sigma: blk_calibration(sigma, data_calib_USA, F0_USA, B_USA, date_settlement),
                     x0, method='Nelder-Mead').x

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

# Assuming you have already defined conversion constants like conv_ACT365

# Convert dates to datetime objects
data_calib_USA_dates = [datetime.strptime(date, '%Y-%m-%d') for date in data_calib_USA['datesExpiry']]
data_calib_EU_dates = [datetime.strptime(date, '%Y-%m-%d') for date in data_calib_EU['datesExpiry']]
date_settlement = datetime.strptime(date_settlement, '%Y-%m-%d')


# Calculate rates for USA and EU
rate_USA, TTM_USA = interp_pricing_params(data_calib_USA_dates, B_USA, date_settlement, year_to_maturity)
rate_EU, TTM_EU = interp_pricing_params(data_calib_EU_dates, B_EU, date_settlement, year_to_maturity)

# Calculate year fractions for USA
yf = np.array([(date - date_settlement).days / conv_ACT365 for date in dates_USA])

# Calculate interpolation date for EU
interp_date = (date_settlement - timedelta(days=1) + timedelta(days=365 * year_to_maturity)).date()

# Extract spot prices
S0_USA = data_USA['spot']
S0_EU = data_EU['spot']

# Compute discount at 1 year
B0_Levy = np.exp(-rate_USA * TTM)

# Concatenate spot prices
S0_Levy = np.array([S0_USA, S0_EU])

# Concatenate rates
rates_Levy = np.array([rate_USA, rate_EU])

# Perform stock simulation
St_Levy = stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0_Levy, rates_Levy, TTM)

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

# BLACK calibration

# Computation of the rates and the time to maturity
rate_USA, TTM = interp_pricing_params(np.array([datenum(date) for date in data_calib_USA['datesExpiry']]), B_USA, date_settlement, year_to_maturity)
rate_EU, _ = interp_pricing_params(np.array([datenum(date) for date in data_calib_EU['datesExpiry']]), B_EU, date_settlement, year_to_maturity)

# Stock prices
S0_USA = data_USA['spot']
S0_EU = data_EU['spot']

# Forward prices
F01_USA = S0_USA * np.exp(rate_USA * TTM)
F01_EU = S0_EU * np.exp(rate_EU * TTM)

# Simulation of the underlying stock prices
St_Black, St_Black_AV = stock_simulation_Black(np.array([sigma_USA, sigma_EU]), np.array([F01_USA, F01_EU]), 
                                                np.array([rate_USA, rate_EU]), rho_historical, TTM)

# Computation of the discount at 1y
B0_black = np.exp(-rate_USA * TTM)


# closed formula

price_semiclosed = blk_semiclosed(data_USA.spot, rate_USA, rate_EU, sigma_USA, sigma_EU, rho_historical, TTM)


#Display of the prices:

disp_contract_prices(mean_price_Levy,IC_Levy,mean_price_Black,IC_Black,mean_price_Black_AV,IC_Black_AV,price_semiclosed)

# End run time computation and display it
end_time = time.time()
print(f"Run time: {end_time - start_time} seconds")
