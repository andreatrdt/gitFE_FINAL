# Library of functions dedicted to the calibration of the model
# 
# USES:
# hist_corr()
# OTM_preprocessing()
# yearfrac()
# blsdelta()
# black_price()
# implied_volatility()
# dataset_preprocessing()
# callPriceLewis_pref()
# RMSE_total()
# calibration()
# nonlinconstr()
# nonlinconstr_corr()
# marginal_param()
# blk_calibration()
#
# Authors :
# M. Maspes
# A. Tarditi
# M. Torba

#%%

# Built in libraries

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


# Requirements for the Black computation
import math
from scipy.stats import norm
from scipy.optimize import brentq, fsolve

# Requirements for the FFT computation
from dateutil.relativedelta import relativedelta
from scipy.integrate import quad
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
from scipy.optimize import OptimizeResult
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

#%%

#####################################################################################################
#################################          CORRELATION FUNCTIONS          ###########################
#####################################################################################################

# Historical correlation
def hist_corr(returns_data):
    returns_EU = returns_data[:,0]
    returns_USA = returns_data[:,1]
    
    rho = np.corrcoef(returns_EU, returns_USA)[0, 1]
    
    return rho

#####################################################################################################
##############################         DATA PROCESSING FUNCTIONS          ###########################
#####################################################################################################


# OTM preprocessing
def OTM_preprocessing(dataset, F0):
    # Introduction
    length_dataset = dataset._length_dates()

    for ii in range(length_dataset):

        # Quantities of interest
        strike_ATM = F0[ii]
        strikes = dataset._strikes_i(ii)
        
        # Compute the related index for the strike ATM
        indicator_call_OTM = strikes > strike_ATM
        indicator_put_OTM = strikes <= strike_ATM

        # Restructuring of the dataset
        dataset._set_callBid(ii, dataset._callBid_i(ii) * indicator_call_OTM)
        dataset._set_callAsk(ii, dataset._callAsk_i(ii) * indicator_call_OTM)
        dataset._set_putBid(ii, dataset._putBid_i(ii) * indicator_put_OTM)
        dataset._set_putAsk(ii, dataset._putAsk_i(ii) * indicator_put_OTM)
    
    return dataset


def yearfrac(start_date, end_date, basis):
    # ACT/365 convention
    if basis == 3:
        return (end_date - start_date).days / 365.0
    else:
        raise NotImplementedError("Other day count conventions are not implemented.")
    

# REFINEMENT DATASET
def blsdelta(S, K, r, T, sigma):
    S = np.array(S, dtype=np.float64)
    K = np.array(K, dtype=np.float64)
    T = np.array(T, dtype=np.float64)
    sigma = np.array(sigma, dtype=np.float64)
    
    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    delta_call = norm.cdf(d1)
    delta_put = delta_call - 1
    return delta_call, delta_put

def black_price(option_type, F, K, T, r, sigma):
    d1 = (math.log(F / K) + (0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)

    if option_type == 'call':
        price = math.exp(-r * T) * (F * norm.cdf(d1) - K* norm.cdf(d2))
    elif option_type == 'put':
        price = math.exp(-r * T) * (K * norm.cdf(-d2) - F * norm.cdf(-d1))
    else:
        raise ValueError("Invalid option type. Use 'call' or 'put'.")

    return price

def implied_volatility(option_type, F, K, T, r, market_price):
    def objective_function(sigma):
        return black_price(option_type, F, K, T, r, sigma) - market_price

    # Brent's method requires an initial interval [a, b] where the function changes sign
    implied_vol = brentq(objective_function, 1e-6, 3.0)  # starting interval for volatility
    return implied_vol

def dataset_preprocessing(dataset, F0, B0, date_settlement):
    # Introduction
    conv_ACT365 = 3
    spot_ATM = dataset._spot
    length_dates = dataset._length_dates()

    for ii in range(length_dates):

        # Set the ATM strike price
        strike_ATM = F0[ii]

        # Extract the strikes from the dataset
        strikes = dataset._strikes_i(ii)
        
        # Calculate time to maturity and interest rate
        TTM = yearfrac(date_settlement, dataset._list_dates_i(ii), conv_ACT365)
        interest_rate = -np.log(B0[ii]) / TTM

        # Calculate the index of the OTM options
        idx_call_OTM = np.where(strikes > strike_ATM)[0]
        idx_put_OTM = np.where(strikes <= strike_ATM)[0]

        # Calculate mid prices for calls and puts
        mid_price_call = (dataset._callAsk_i(ii) + dataset._callBid_i(ii)) / 2
        mid_price_call = mid_price_call[idx_call_OTM]

        mid_price_put = (dataset._putAsk_i(ii) + dataset._putBid_i(ii)) / 2
        mid_price_put = mid_price_put[idx_put_OTM]

        # Computation of the implied volatilities
        impvol_call_i = np.zeros(len(idx_call_OTM))
        impvol_put_i = np.zeros(len(idx_put_OTM))

        for jj in idx_call_OTM:
            impvol_call_i[jj - len(idx_put_OTM)] = (
                implied_volatility('call', F0[ii], strikes[jj], TTM, interest_rate, mid_price_call[jj - len(idx_put_OTM)]))

        for jj in idx_put_OTM:
            impvol_put_i[jj] = (
                implied_volatility('put', F0[ii], strikes[jj], TTM, interest_rate, mid_price_put[jj]))

        # Computation of the deltas
        delta_call, _ = blsdelta(spot_ATM, strikes[idx_call_OTM], interest_rate, TTM, impvol_call_i)
        _, delta_put = blsdelta(spot_ATM, strikes[idx_put_OTM], interest_rate, TTM, impvol_put_i)

        indicator_put = (delta_put <= -0.1) & (delta_put >= -0.9)
        idx_put = np.where(indicator_put)[0]
        
        indicator_call = (delta_call >= 0.1) & (delta_call <= 0.9)
        idx_call = np.where(indicator_call)[0] + len(delta_put)

        # Restructuring of the dataset
        dataset._set_callBid(ii, dataset._callBid_i(ii)[idx_call])
        dataset._set_callAsk(ii, dataset._callAsk_i(ii)[idx_call])
        dataset._set_putBid(ii, dataset._putBid_i(ii)[idx_put])
        dataset._set_putAsk(ii, dataset._putAsk_i(ii)[idx_put])

        idx_combined = np.unique(np.concatenate((idx_put, idx_call)))
        dataset._set_strikes(ii, dataset._strikes_i(ii)[idx_combined])
    
    return dataset


#####################################################################################################
#################################          CALIBRATION FUNCTIONS          ############################
#####################################################################################################

def callPriceLewis_pref(B0, F0, log_moneyness, sigma, k, theta, TTM, M, dz):
    # Introduction of the FFT parameters
    
    # Initial constraints for FFT parameters
    N = 2 ** M

    dx = (2 * np.pi) / (N * dz)
    xN = ((N - 1) * dx) / 2
    x1 = -xN
    zN = ((N - 1) * dz) / 2
    z1 = -zN

    # Computation of the grid
    X = np.linspace(x1, xN, N)
    Z = np.linspace(z1, zN, N)

    # Computation of prefactor
    preFactor = dx * np.exp(-1j * x1 * Z)

    # Computation of FFT input
    xi = -X - 1j / 2

    # Computation of the characteristic function
    charFct = np.exp(TTM * (1 / k * (1 - np.sqrt(1 - 2j * xi * k * theta + xi ** 2 * k * sigma ** 2))) - \
                    xi * 1j * TTM * 1 / k * (1 - np.sqrt(1 - 2 * k * theta - k * sigma ** 2)))

    # Computation of the integrand
    integrand = 1 / (2 * np.pi) * 1 / (X ** 2 + 1 / 4) * charFct

    j = np.arange(0, N)
    inputFFT = integrand * np.exp(-1j * z1 * dx * j)

    # Computation of the integral
    IntLewis = preFactor * fft(inputFFT)

    IntLewis = np.real(IntLewis)
    IntLewis = interp1d(Z, IntLewis, kind='linear', fill_value='extrapolate')(log_moneyness)

    # Computation of the Call Price through Lewis
    price = B0 * F0 * (1 - np.exp(-log_moneyness / 2) * IntLewis)

    return price

# compute the RMSE
def RMSE_total(params, dataset, F0, B0, date_settlement):

    # Unpack the parameters
    k, theta, sigma = params

    # Conventions
    conv_ACT365 = 3

    # Initial parameters
    RMSE = np.zeros(dataset._length_dates())
    N_options = 0
    
    # FFT parameters
    M = 15
    dz = 0.001

    # Computation 
    for ii in range(dataset._length_dates()):

        # Initialization
        put_length = len(dataset._putBid_i(ii))

        N_options += len(dataset._strikes_i(ii))

        # Logmoneyness values
        log_moneyness = np.log(F0[ii] / dataset._strikes_i(ii))

        # Time to maturity
        TTM = yearfrac(date_settlement, dataset._list_dates_i(ii), conv_ACT365)

        # Pricing 
        prices = callPriceLewis_pref(B0[ii], F0[ii], log_moneyness, sigma, k, theta, TTM, M, dz)

        call_prices = prices[put_length:]
        put_prices = prices[:put_length] - F0[ii] * B0[ii] + np.array(dataset._strikes_i(ii))[:put_length] * B0[ii]

        mean_call_price = (dataset._callAsk_i(ii) + dataset._callBid_i(ii)) / 2
        mean_put_price = (dataset._putAsk_i(ii) + dataset._putBid_i(ii)) / 2

        # Computation of RMSE
        RMSE[ii] = np.sqrt(np.sum((call_prices - mean_call_price) ** 2) + np.sum((put_prices - mean_put_price) ** 2))

    # Final adjusting of RMSE
    RMSE_total = np.sum(RMSE)

    return RMSE_total

# Calibration function
def calibration(params, data_EU, data_USA, F0_EU, B0_EU, F0_USA, B0_USA, date_settlement):
    # Unpacking of the parameters
    params_USA = params[:3]
    params_EU = params[3:]

    # Computation of the RMSE
    RMSE_EU = RMSE_total(params_EU, data_EU, F0_EU, B0_EU, date_settlement)
    RMSE_USA = RMSE_total(params_USA, data_USA, F0_USA, B0_USA, date_settlement)

    # Computation of the weights
    weight_USA = data_USA._spot/ (data_EU._spot + data_USA._spot)
    weight_EU = data_EU._spot/ (data_EU._spot + data_USA._spot)

    # Computation of the distance
    distance = weight_EU * RMSE_EU + weight_USA * RMSE_USA

    return distance

# non linear constraints
def nonlinconstr(x):
    # Unpacking the parameters
    k1, theta1, sigma1, k2, theta2, sigma2 = x

    # Constraints on the equalities
    # Constraint to create the entire equality given on the final parameters
    ceq = (sigma1**2 / (k1 * theta1**2)) - (sigma2**2 / (k2 * theta2**2))

    # Constraints on the inequalities
    c = []

    return c, ceq

# non linear constraints for the correlation
def nonlinconstr_corr(params, k1, k2):

    # Unpacking the constraints
    nu_1, nu_2, nu_z = params

    # Constraints on the equalities
    ceq = [
        nu_1 * nu_z / (nu_1 + nu_z) - k1,
        nu_2 * nu_z / (nu_2 + nu_z) - k2]


    # (nu_1 * nu_2 / ((nu_1 + nu_z) * (nu_2 + nu_z))) ** 0.5 - ((k1 * k2) ** 0.5 / nu_z)

    # Constraints on the inequalities
    c = []

    return c, ceq

# marginal parameters
def marginal_param(params_USA, params_EU, nu_z):

    # sole system of equations
    def equations(x):
        kappa_1, theta_1, sigma_1 = params_USA
        kappa_2, theta_2, sigma_2 = params_EU
        a_1, a_2, Beta_z, gamma_z = x

        eq1 = a_1 * Beta_z - (kappa_1 * theta_1 / nu_z)
        eq2 = a_2 * Beta_z - (kappa_2 * theta_2 / nu_z)
        eq3 = kappa_1 * sigma_1**2 - nu_z * a_1**2 * gamma_z**2
        eq4 = kappa_2 * sigma_2**2 - nu_z * a_2**2 * gamma_z**2

        return [eq1, eq2, eq3, eq4]

    x0 = 0.2 * np.ones(4)
    sol = fsolve(equations, x0)
    
    return sol

#####################################################################################################
#################################          BLACK FUNCTIONS          ###############################
#####################################################################################################

def blk_calibration(sigma, dataset, F0, B0, date_settlement):

    # Conventions
    conv_ACT365 = 3

    # Initial parameters
    RMSE = np.zeros(dataset._length_dates())

    # Computation
    for ii in range(dataset._length_dates()):
        # Initialization
        put_length = len(dataset._putAsk_i(ii))
        strikes = dataset._strikes_i(ii)

        # Time to maturity
        TTM = yearfrac(date_settlement, dataset._list_dates_i(ii), conv_ACT365)
        
        # Interest rate
        interest_rate = -np.log(B0[ii]) / TTM

        prices = np.zeros(len(strikes))

        # Pricing
        # Price of Call/Puts through the Black formula
        for jj in range(len(strikes)):
            prices[jj] = black_price('call', F0[ii], strikes[jj], TTM, interest_rate, sigma)
        
        call_prices = prices[put_length:]
        put_prices = prices[:put_length] - F0[ii] * B0[ii] + strikes[:put_length] * B0[ii]
        
        mean_call_price = (dataset._callAsk_i(ii) + dataset._callBid_i(ii)) / 2
        mean_put_price = (dataset._putAsk_i(ii) + dataset._putBid_i(ii)) / 2

        # Computation of RMSE
        RMSE[ii] = np.sqrt(np.mean((np.concatenate([call_prices, put_prices]) - np.concatenate([mean_call_price, mean_put_price])) ** 2))

    # Final adjusting of RMSE
    return np.sum(RMSE)