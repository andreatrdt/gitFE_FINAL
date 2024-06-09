# Library of function used forthe ricing of our contract : max( St_USA - S0_USA, 0 ) * I{St_EU < (0.95 * S0_EU)}
#
# USES:
# interp_pricing_params()
# rate_interpolation()
# interp_pricing_params()
# stock_simulation_Levy()
# stock_simulation_Black()
#
# Authors :
# M. Maspes
# A. Tarditi
# M. Torba

#%% 

# Built in libraries
# import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.interpolate import interp1d


from scipy.stats import norm
from dateutil.relativedelta import relativedelta
from scipy.integrate import quad
from scipy.fftpack import fft, ifft

from scipy.interpolate import interp1d
from scipy.optimize import OptimizeResult

from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.stats import multivariate_normal

import calibration_functions as calib_fun
import forward_functions as fwd_fun

# Function to interpolate interest rates
def interp_pricing_params(dates, B0, date_settlement, year_to_maturity):
    conv_ACT365 = 3  # Assuming ACT/365 convention

    # Calculate interpolation date
    interp_date = date_settlement + timedelta(days=365 * year_to_maturity)
    print(f"Interpolation Date: {interp_date}")

    # Calculate time to maturity
    TTM = calib_fun.yearfrac(date_settlement, interp_date, conv_ACT365)
    print(f"Time to Maturity: {TTM}")

    # Calculate interest rate
    rate = rate_interpolation(dates, B0, date_settlement, interp_date)
    print(f"Interest Rate: {rate}")

    return rate, TTM

# Function to interpolate interest rates
def rate_interpolation(dates, B0, settlement_date, interp_date):

    yf = np.zeros(len(dates))

    for ii in range(len(dates)):
        yf[ii] = calib_fun.yearfrac(settlement_date, dates[ii], 3)

    TTM = calib_fun.yearfrac(settlement_date, interp_date, 3)

    rates = -np.log(B0) / yf
    rate_interpolator = interp1d(yf, rates)

    return rate_interpolator(TTM)

# Function to extract interest rate and time to maturity
def interp_pricing_params(dates, B0, date_settlement, year_to_maturity):
    conv_ACT365 = 3  # Assuming ACT/365 convention

    # Calculate interpolation date
    interp_date = (date_settlement  + timedelta(days=366 * year_to_maturity))

    # Calculate time to maturity
    TTM = calib_fun.yearfrac(date_settlement, interp_date, conv_ACT365)

    # Calculate interest rate
    rate = rate_interpolation(dates, B0, date_settlement, interp_date)

    return rate, TTM

# Function to simulate stock prices
def stock_simulation_Levy(idiosync_USA, idiosync_EU, syst_Z, params_USA, params_EU, S0, rates, TTM):

    # Unpacking parameters
    kappa_USA, theta_USA, sigma_USA = params_USA
    kappa_EU, theta_EU, sigma_EU = params_EU

    nu_z = syst_Z[0]
    nu_USA, nu_EU = idiosync_USA[0], idiosync_EU[0]

    Beta_z = syst_Z[1]
    Beta_USA, Beta_EU = idiosync_USA[1], idiosync_EU[1]

    gamma_z = syst_Z[2]
    gamma_USA, gamma_EU = idiosync_USA[2], idiosync_EU[2]

    a_USA, a_EU = idiosync_USA[3], idiosync_EU[3]
    
    # Computation of the support parameters
    nSim = 10000000
    
    drift_compensator_USA = - 1/kappa_USA * (1 - np.sqrt(1 - 2*kappa_USA*theta_USA - kappa_USA*sigma_USA**2))
    drift_compensator_EU = - 1/kappa_EU * (1 - np.sqrt(1 - 2*kappa_EU*theta_EU - kappa_EU*sigma_EU**2))
    drift_compensator = np.array([drift_compensator_USA, drift_compensator_EU])
    
    # Simulation of the NIG process
    g_1 = np.random.randn(nSim)
    g_2 = np.random.randn(nSim)
    g_z = np.random.randn(nSim)
    
    G_1 = np.random.wald(1, TTM/nu_USA, size=nSim)
    G_2 = np.random.wald(1, TTM/nu_EU, size=nSim)
    G_z = np.random.wald(1, TTM/nu_z, size=nSim)
    
    # Y_1 = -(0.5 + Beta_USA) * gamma_USA**2 * G_1 * TTM + gamma_USA * np.sqrt(TTM * G_1) * g_1
    # Y_2 = -(0.5 + Beta_EU) * gamma_EU**2 * G_2 * TTM + gamma_EU * np.sqrt(TTM * G_2) * g_2
    # Z = -(0.5 + Beta_z) * gamma_z**2 * G_z * TTM + gamma_z * np.sqrt(TTM * G_z) * g_z

    Y_1 = Beta_USA * G_1 * TTM + gamma_USA * np.sqrt(TTM * G_1) * g_1
    Y_2 = Beta_EU * G_2 * TTM + gamma_EU * np.sqrt(TTM * G_2) * g_2
    Z = Beta_z * G_z * TTM + gamma_z * np.sqrt(TTM * G_z) * g_z

    # Marginal processes
    X_1 = Y_1 + a_USA * Z
    X_2 = Y_2 + a_EU * Z
    
    Xt = np.column_stack((X_1, X_2))
    
    stock = S0 * np.exp((rates + drift_compensator) * TTM + Xt)
    
    return stock

# Function to simulate stock prices
def stock_simulation_Black(sigmas, F0, rates, rho, TTM):
    nSim = 10000000

    # Simulation of the NIG process
    meanVector = np.array([0, 0])
    covarianceMatrix = np.array([[1, rho * 1], [rho * 1, 1]])

    g = multivariate_normal.rvs(mean=meanVector, cov=covarianceMatrix, size=nSim)

    # Creation of Xt dynamic
    Xt = (-0.5 * sigmas ** 2 * TTM + sigmas * np.sqrt(TTM) * g.T).T

    # Computation of the initial stock
    prices = F0 * np.exp(Xt)

    return prices