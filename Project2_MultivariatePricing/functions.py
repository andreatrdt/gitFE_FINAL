import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import scipy.io as sio
import datetime
import warnings
import os
from scipy.stats import norm
from dateutil.relativedelta import relativedelta
import numpy as np
from scipy.integrate import quad
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
from scipy.optimize import OptimizeResult
from datetime import datetime, timedelta
from scipy.interpolate import interp1d




# dataset exploration
def dataset_exploration(data_EU, data_USA, date_settlement):
    """
    Exploration of the dataset and consequent plots
    
    INPUT:
    data_EU:               dataset EU
    data_USA:              dataset USA
    date_settlement:       initial date of the computation
    """

    # Plot of the European Call Prices
    plt.figure()
    for ii in range(len(data_EU['datesExpiry'])):
        # Compute the mid prices
        mid_call = (data_EU['callAsk'][ii]['prices'] + data_EU['callBid'][ii]['prices']) / 2
        plt.plot(data_EU['strikes'][ii]['value'], mid_call, label=f'Expiry {ii+1}')
    plt.grid(True)
    plt.title('European Call Prices')
    plt.xlabel('Strikes')
    plt.ylabel('Prices')
    plt.legend()
    plt.show()

    # Plot of the European Put Prices
    plt.figure()
    for ii in range(len(data_EU['datesExpiry'])):
        # Compute the mid prices
        mid_put = (data_EU['putAsk'][ii]['prices'] + data_EU['putBid'][ii]['prices']) / 2
        plt.plot(data_EU['strikes'][ii]['value'], mid_put, label=f'Expiry {ii+1}')
    plt.grid(True)
    plt.title('European Put Prices')
    plt.xlabel('Strikes')
    plt.ylabel('Prices')
    plt.legend()
    plt.show()

    # Plot of the American Call Prices
    plt.figure()
    for ii in range(len(data_USA['datesExpiry'])):
        # Compute the mid prices
        mid_call = (data_USA['callAsk'][ii]['prices'] + data_USA['callBid'][ii]['prices']) / 2
        plt.plot(data_USA['strikes'][ii]['value'], mid_call, label=f'Expiry {ii+1}')
    plt.grid(True)
    plt.title('American Call Prices')
    plt.xlabel('Strikes')
    plt.ylabel('Prices')
    plt.legend()
    plt.show()

    # Plot of the American Put Prices
    plt.figure()
    for ii in range(len(data_USA['datesExpiry'])):
        # Compute the mid prices
        mid_put = (data_USA['putAsk'][ii]['prices'] + data_USA['putBid'][ii]['prices']) / 2
        plt.plot(data_USA['strikes'][ii]['value'], mid_put, label=f'Expiry {ii+1}')
    plt.grid(True)
    plt.title('American Put Prices')
    plt.xlabel('Strikes')
    plt.ylabel('Prices')
    plt.legend()
    plt.show()

    # Study of the penny options
    # Penny options are those options with value lower than 0.1 index points

    # European
    count_EU_call_ask = sum([np.sum(prices <= 0.1) for prices in data_EU['callAsk'].flatten()])
    count_EU_call_bid = sum([np.sum(prices <= 0.1) for prices in data_EU['callBid'].flatten()])
    count_EU_put_ask = sum([np.sum(prices <= 0.1) for prices in data_EU['putAsk'].flatten()])
    count_EU_put_bid = sum([np.sum(prices <= 0.1) for prices in data_EU['putBid'].flatten()])

    print('The number of European Penny options is')
    print(count_EU_call_ask + count_EU_call_bid + count_EU_put_ask + count_EU_put_bid)

    # American
    count_USA_call_ask = sum([np.sum(prices <= 0.1) for prices in data_USA['callAsk'].flatten()])
    count_USA_call_bid = sum([np.sum(prices <= 0.1) for prices in data_USA['callBid'].flatten()])
    count_USA_put_ask = sum([np.sum(prices <= 0.1) for prices in data_USA['putAsk'].flatten()])
    count_USA_put_bid = sum([np.sum(prices <= 0.1) for prices in data_USA['putBid'].flatten()])

    print('The number of American Penny options is')
    print(count_USA_call_ask + count_USA_call_bid + count_USA_put_ask + count_USA_put_bid)

    # Study of the liquidity criterion

    # European
    indicator_liq_call = 0
    indicator_liq_put = 0
    for ii in range(len(data_EU['datesExpiry'])):
        liquidity_call = (data_EU['callAsk'][ii]['prices'] - data_EU['callBid'][ii]['prices']) / data_EU['callAsk'][ii]['prices']
        indicator_liq_call += np.sum(liquidity_call >= 0.6)

        liquidity_put = (data_EU['putAsk'][ii]['prices'] - data_EU['putBid'][ii]['prices']) / data_EU['putAsk'][ii]['prices']
        indicator_liq_put += np.sum(liquidity_put >= 0.6)

    print('The number of illiquid options in the EU market is:')
    print(indicator_liq_call + indicator_liq_put)

    # American
    indicator_liq_call = 0
    indicator_liq_put = 0
    for ii in range(len(data_USA['datesExpiry'])):
        liquidity_call = (data_USA['callAsk'][ii]['prices'] - data_USA['callBid'][ii]['prices']) / data_USA['callAsk'][ii]['prices']
        indicator_liq_call += np.sum(liquidity_call >= 0.6)

        liquidity_put = (data_USA['putAsk'][ii]['prices'] - data_USA['putBid'][ii]['prices']) / data_USA['putAsk'][ii]['prices']
        indicator_liq_put += np.sum(liquidity_put >= 0.6)

    print('The number of illiquid options in the USA market is:')
    print(indicator_liq_call + indicator_liq_put)
    
# Historical correlation
def hist_corr(dataset):
    # Extract returns data for USA and EU from the dataset
    returns_USA = dataset[:, 0]
    returns_EU = dataset[:, 1]
    
    # Compute the correlation coefficient
    rho = np.corrcoef(returns_EU, returns_USA)[0, 1]
    
    return rho

# discount factor
def estimation_discount_factor(Gi, Ki):
    """
    Computation of the estimated discount factor B_bar(0, T)
    
    INPUT:
    Gi:               [ARRAY] synthetic forwards for the estimation
    Ki:               [ARRAY] strikes for the estimation
    
    OUTPUT:
    B_bar:            estimated discount factor for the forward computation
    """
    
    # Computation of the mean values
    # Mean values of the strikes and synthetic forwards
    G_hat = np.mean(Gi)
    K_hat = np.mean(Ki)
    
    # Final computation of B_bar(0, T)
    B_bar = -np.dot(Ki - K_hat, Gi - G_hat) / np.dot(Ki - K_hat, Ki - K_hat)
    
    return B_bar


#syntehtic forward
def synthetic_forward(call_bid, call_ask, put_bid, put_ask):
    """
    Computation of the Synthetic Forward Gi for each of the required horizons,
    the SF is the building block of the Forward computation

    INPUT:
    call_bid:          [ARRAY] bid prices of the Call
    call_ask:          [ARRAY] ask prices of the Call
    put_bid:           [ARRAY] bid prices of the Put
    put_ask:           [ARRAY] ask prices of the Put

    OUTPUT:
    Gi:                Synthetic Forward for each strike
    Gi_ask:            Synthetic Forward Ask for each strike
    Gi_bid:            Synthetic Forward Bid for each strike
    """
    
    # Computation Synthetic Forward Bid
    Gi_bid = call_bid - put_ask

    # Computation Synthetic Forward Ask
    Gi_ask = call_ask - put_bid

    # Computation final Synthetic Forward
    Gi = (Gi_bid + Gi_ask) / 2

    return Gi, Gi_ask, Gi_bid


# forward prices 
def forward_prices(dataset, flag):
    """
    Computation of the forward prices following the Baviera, Azzone paper
    
    INPUT:
    dataset:           [DICT] data containing all the required tables
    flag:              [0: without plots & slope; 1: with plots & slope]
    
    OUTPUT:
    F_vector:          [ARRAY] forwards value F(0, T)
    B_bar_vector:      [ARRAY] market calibrated discount B(0, T)
    
    USES:
    function synthetic_forward()
    function estimation_discount_factor()
    """

    # Introduction of the return vectors
    F_vector = np.zeros(len(dataset['datesExpiry']))
    B_bar_vector = np.zeros(len(dataset['datesExpiry']))

    # Computation of the forwards
    for ii in range(len(dataset['datesExpiry'])):

        # Computation of the strikes and synthetic forwards
        Ki = dataset['strikes'][ii]['value']

        # Synthetic forwards
        Gi, Gi_ask, Gi_bid = synthetic_forward(dataset['callBid'][ii]['prices'], dataset['callAsk'][ii]['prices'],
                                               dataset['putBid'][ii]['prices'], dataset['putAsk'][ii]['prices'])

        # Computation of the estimated discount factor
        B_bar_vector[ii] = estimation_discount_factor(Gi, Ki)

        # Computation of the forward prices
        F_i_vector = Gi / B_bar_vector[ii] + Ki
        F_ask__i_vector = Gi_ask / B_bar_vector[ii] + Ki
        F_bid__i_vector = Gi_bid / B_bar_vector[ii] + Ki

        # Computation of the required Forwards
        F_vector[ii] = np.mean(F_i_vector)

        # Eventual plot
        if flag:
            plt.figure()
            plt.plot(Ki, F_ask__i_vector, '*-', label='Ask')
            plt.plot(Ki, F_i_vector, 'o-', label='Mid')
            plt.plot(Ki, F_bid__i_vector, '*-', label='Bid')
            plt.grid(True)
            plt.xlabel('Strikes')
            plt.ylabel('Forward Prices')
            plt.legend()
            plt.title(f'Forward Prices at Expiry {ii+1}')
            plt.show()

            # Computation of the slope
            p = np.polyfit(Ki, F_i_vector - Ki, 1)
            slope = p[0]
            print(f'Estimated Slope for Expiry {ii+1}: {slope}')

    return F_vector, B_bar_vector


# OTM preprocessing
def OTM_preprocessing(dataset, F0):
    """
    Computation of the OTM options only
    
    INPUT:
    dataset:          [DICT] initial dataset
    F0:               [ARRAY] forward price at time 0
    
    OUTPUT:
    dataset:          [DICT] dataset modified for the OTM options 
    """
    for ii in range(len(dataset['datesExpiry'])):
        # Quantities of interest
        strike_ATM = F0[ii]
        strikes = dataset['strikes'][ii]['value']
        
        # Compute the related index for the strike ATM
        indicator_call_OTM = strikes > strike_ATM
        indicator_put_OTM = strikes <= strike_ATM
        
        # Restructuring of the dataset
        dataset['callBid'][ii]['prices'] *= indicator_call_OTM
        dataset['callAsk'][ii]['prices'] *= indicator_call_OTM
        dataset['putBid'][ii]['prices'] *= indicator_put_OTM
        dataset['putAsk'][ii]['prices'] *= indicator_put_OTM
    
    return dataset


def yearfrac(start_date, end_date, basis):
    # ACT/365 convention
    if basis == 3:
        return (end_date - start_date).days / 365.0
    else:
        raise NotImplementedError("Other day count conventions are not implemented.")

def blkimpv(F, K, r, T, market_price, opt_type='Call'):
    """
    Compute Black implied volatility. This is a placeholder function.
    The actual implementation should use a numerical solver.
    """
    # This is a simplified placeholder for demonstration purposes.
    # In practice, you'd use a numerical method to solve for implied volatility.
    sigma = np.full_like(K, 0.2)  # Placeholder constant volatility
    return sigma

def blsdelta(S, K, r, T, sigma):
    """
    Compute Black-Scholes delta. Placeholder function using simplified calculations.
    """
    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    delta_call = norm.cdf(d1)
    delta_put = delta_call - 1
    return delta_call, delta_put

def dataset_preprocessing(dataset, F0, B0, date_settlement, flag):
    """
    Preprocessing of the dataset for each maturity in order to remove all
    those options without delta_Black between [10%, 90%]
    
    INPUT:
    dataset:            [DICT] initial dataset
    F0:                 [ARRAY] initial forward F(0, T)
    B0:                 [ARRAY] initial discount B(0, T)
    date_settlement:    [DATETIME] settlement date
    flag:               [0: no plots, 1: with plots]
    
    OUTPUT:
    dataset:            modified version of the dataset
    """
    conv_ACT365 = 3

    for ii in range(len(dataset['datesExpiry'])):
        spot_ATM = dataset['spot']
        strike_ATM = F0[ii]
        strikes = dataset['strikes'][ii]['value']
        
        TTM = yearfrac(date_settlement, dataset['datesExpiry'][ii], conv_ACT365)
        interest_rate = -np.log(B0[ii]) / TTM

        idx_call_OTM = np.where(strikes > strike_ATM)[0]
        idx_put_OTM = np.where(strikes <= strike_ATM)[0]

        mid_price_call = (dataset['callAsk'][ii]['prices'] + dataset['callBid'][ii]['prices']) / 2
        mid_price_call = mid_price_call[idx_call_OTM]
        
        mid_price_put = (dataset['putAsk'][ii]['prices'] + dataset['putBid'][ii]['prices']) / 2
        mid_price_put = mid_price_put[idx_put_OTM]

        if len(dataset['datesExpiry']) == 20 and ii == len(dataset['datesExpiry']) - 1:
            strikes_removal_put = [2000, 2200, 2400, 2600, 2700, 4900, 5000, 5100, 5200]
            strikes_removal_call = [5300, 7200]
            
            idx_removed_strikes_put = np.isin(strikes[idx_put_OTM], strikes_removal_put)
            idx_removed_strikes_call = np.isin(strikes[idx_call_OTM], strikes_removal_call)
            
            idx_put = np.where(~idx_removed_strikes_put)[0]
            idx_call = np.where(~idx_removed_strikes_call)[0] + len(idx_put_OTM)
            
            strikes = np.concatenate((strikes[idx_put], strikes[idx_call]))
            
            mid_price_put = mid_price_put[idx_put]
            mid_price_call = mid_price_call[idx_call - len(idx_put_OTM)]
            
            idx_call_OTM = np.where(strikes > strike_ATM)[0]
            idx_put_OTM = np.where(strikes <= strike_ATM)[0]
            
            mid_price_call = mid_price_call[idx_call_OTM - len(idx_put_OTM)]
            mid_price_put = mid_price_put[idx_put_OTM]

        impvol_call_i = blkimpv(F0[ii], strikes[idx_call_OTM], interest_rate, TTM, mid_price_call, 'Call')
        impvol_put_i = blkimpv(F0[ii], strikes[idx_put_OTM], interest_rate, TTM, mid_price_put, 'Put')

        if flag:
            plt.figure()
            plt.plot(strikes[idx_call_OTM], impvol_call_i, 'o-', label='Implied vol Call')
            plt.plot(strikes[idx_put_OTM], impvol_put_i, '*-', label='Implied vol Put')
            plt.grid(True)
            plt.title('Implied volatilities')
            plt.xlabel('Strikes')
            plt.ylabel('Volatilities')
            plt.legend()
            plt.show()

        delta_call, _ = blsdelta(spot_ATM, strikes[idx_call_OTM], interest_rate, TTM, impvol_call_i)
        _, delta_put = blsdelta(spot_ATM, strikes[idx_put_OTM], interest_rate, TTM, impvol_put_i)

        indicator_put = (delta_put <= -0.1) & (delta_put >= -0.9)
        idx_put = np.where(indicator_put)[0]
        
        indicator_call = (delta_call >= 0.1) & (delta_call <= 0.9)
        idx_call = np.where(indicator_call)[0] + len(delta_put)
        
        dataset['callBid'][ii]['prices'] = dataset['callBid'][ii]['prices'][idx_call]
        dataset['callAsk'][ii]['prices'] = dataset['callAsk'][ii]['prices'][idx_call]
        dataset['callAsk'][ii]['impvol'] = impvol_call_i[idx_call - len(delta_put)]
        
        dataset['putAsk'][ii]['prices'] = dataset['putAsk'][ii]['prices'][idx_put]
        dataset['putBid'][ii]['prices'] = dataset['putBid'][ii]['prices'][idx_put]
        dataset['putBid'][ii]['impvol'] = impvol_put_i[idx_put]

        idx_combined = np.unique(np.concatenate((idx_put, idx_call)))
        dataset['strikes'][ii]['value'] = dataset['strikes'][ii]['value'][idx_combined]
        
        dataset['Volume_call'][ii]['volume'] = dataset['Volume_call'][ii]['volume'][idx_call]
        dataset['Volume_put'][ii]['volume'] = dataset['Volume_put'][ii]['volume'][idx_put]
    
    return dataset


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


def RMSE_total(params, dataset, F0, B0, date_settlement):
    # Unpack the parameters
    k, theta, sigma = params

    # Conventions
    conv_ACT365 = 3

    # Initial parameters
    RMSE = np.zeros(len(dataset['datesExpiry']))
    weights = np.zeros(len(dataset['datesExpiry']))
    N_options = 0
    
    # FFT parameters
    M = 15
    dz = 0.001

    # Computation 
    for ii in range(min(len(dataset['datesExpiry']), 19)):
        # Initialization
        put_length = len(dataset['putAsk'][ii]['prices'])

        N_options += len(dataset['strikes'][ii]['value'])

        # Logmoneyness values
        log_moneyness = np.log(F0[ii] / np.array(dataset['strikes'][ii]['value']))

        # Time to maturity
        TTM = (date_settlement - dataset['datesExpiry'][ii]).days / 365.0

        # Pricing 
        prices = callPriceLewis_pref(B0[ii], F0[ii], log_moneyness, sigma, k, theta, TTM, M, dz)

        call_prices = prices[put_length:]
        put_prices = prices[:put_length] - F0[ii] * B0[ii] + np.array(dataset['strikes'][ii]['value'])[:put_length] * B0[ii]

        mean_call_price = (np.array(dataset['callAsk'][ii]['prices']) + np.array(dataset['callBid'][ii]['prices'])) / 2
        mean_put_price = (np.array(dataset['putAsk'][ii]['prices']) + np.array(dataset['putBid'][ii]['prices'])) / 2

        # Computation of RMSE
        RMSE[ii] = np.sqrt(np.sum((call_prices - mean_call_price) ** 2) + np.sum((put_prices - mean_put_price) ** 2))

    # Final adjusting of RMSE
    RMSE_total = np.sum(RMSE)

    return RMSE_total


def new_calibration(params, data_EU, data_USA, F0_EU, B0_EU, F0_USA, B0_USA, date_settlement):
    # Unpacking of the parameters
    params_USA = params[:3]
    params_EU = params[3:]

    # Computation of the RMSE
    RMSE_EU = RMSE_total(params_EU, data_EU, F0_EU, B0_EU, date_settlement)
    RMSE_USA = RMSE_total(params_USA, data_USA, F0_USA, B0_USA, date_settlement)

    # Computation of the weights
    weight_USA = data_USA['spot'] / (data_EU['spot'] + data_USA['spot'])
    weight_EU = data_EU['spot'] / (data_EU['spot'] + data_USA['spot'])

    # Computation of the distance
    distance = weight_EU * RMSE_EU + weight_USA * RMSE_USA

    return distance

def nonlinconstr(x):
    # Unpacking the parameters
    k1, theta1, sigma1, k2, theta2, sigma2 = x

    # Constraints on the equalities
    # Constraint to create the entire equality given on the final parameters
    ceq = (sigma1**2 / (k1 * theta1**2)) - (sigma2**2 / (k2 * theta2**2))

    # Constraints on the inequalities
    c = []

    return c, ceq

def nonlinconstr_corr(params, k1, k2):
    # Unpacking the constraints
    nu_1, nu_2, nu_z = params

    # Constraints on the equalities
    ceq = [
        nu_1 * nu_z / (nu_1 + nu_z) - k1,
        nu_2 * nu_z / (nu_2 + nu_z) - k2,
        (nu_1 * nu_2 / ((nu_1 + nu_z) * (nu_2 + nu_z)))**0.5 - (k1 * k2)**0.5 / nu_z
    ]

    # Constraints on the inequalities
    c = []

    return c, ceq

def marginal_param(params_USA, params_EU, nu_z, rho):
    def equations(x):
        kappa_1, theta_1, sigma_1 = params_USA
        kappa_2, theta_2, sigma_2 = params_EU
        a_1, a_2, Beta_z, gamma_z = x

        eq1 = a_1 * Beta_z - (kappa_1 * theta_1 / nu_z)
        eq2 = a_2 * Beta_z - (kappa_2 * theta_2 / nu_z)
        eq3 = kappa_1 * sigma_1**2 - nu_z * a_1**2 * gamma_z**2
        eq4 = kappa_2 * sigma_2**2 - nu_z * a_2**2 * gamma_z**2

        return [eq1, eq2, eq3, eq4]

    x0 = np.ones(4)
    sol = OptimizeResult()
    sol.x = fsolve(equations, x0)
    
    return sol


def disp_marginal_params(sol_USA, sol_EU, Beta_z, gamma_z, nu_z, save_results):
    # Extract parameters for the USA
    a_USA = sol_USA.x[3]
    Beta_USA = sol_USA.x[1]
    gamma_USA = sol_USA.x[2]
    nu_1 = sol_USA.x[0]

    # Extract parameters for the EU
    a_EU = sol_EU.x[3]
    Beta_EU = sol_EU.x[1]
    gamma_EU = sol_EU.x[2]
    nu_2 = sol_EU.x[0]

    # Display the parameters
    print('PARAMETERS obtained by convolution:')
    print('----------------------------')
    print('a_USA:', a_USA)
    print('a_EU:', a_EU)
    print('----------------------------')
    
    # Display Beta parameters
    print('Beta_z:', Beta_z)
    print('Beta_USA:', Beta_USA)
    print('Beta_EU:', Beta_EU)
    print('----------------------------')
    
    # Display gamma parameters
    print('gamma_z:', gamma_z)
    print('gamma_USA:', gamma_USA)
    print('gamma_EU:', gamma_EU)
    print('----------------------------')

    print('nu_1:', nu_1)
    print('nu_2:', nu_2)
    print('nu_z:', nu_z)
    print('----------------------------')

    if save_results:
        # Save the parameters to a text file
        with open('results.txt', 'a') as file:
            file.write('PARAMETERS obtained by convolution:\n')
            file.write('----------------------------\n')
            file.write('a_USA: {}\n'.format(a_USA))
            file.write('a_EU: {}\n'.format(a_EU))
            file.write('----------------------------\n')
            file.write('Beta_z: {}\n'.format(Beta_z))
            file.write('Beta_USA: {}\n'.format(Beta_USA))
            file.write('Beta_EU: {}\n'.format(Beta_EU))
            file.write('----------------------------\n')
            file.write('gamma_z: {}\n'.format(gamma_z))
            file.write('gamma_USA: {}\n'.format(gamma_USA))
            file.write('gamma_EU: {}\n'.format(gamma_EU))
            file.write('----------------------------\n')
            file.write('nu_1: {}\n'.format(nu_1))
            file.write('nu_2: {}\n'.format(nu_2))
            file.write('nu_z: {}\n'.format(nu_z))
            file.write('----------------------------\n')

def disp_params(params_marginals, initial_cond, flag):
    # Unpack the calibrated parameters
    params_USA = params_marginals[:3]
    params_EU = params_marginals[3:]

    kappa_USA, theta_USA, sigma_USA = params_USA
    kappa_EU, theta_EU, sigma_EU = params_EU

    if flag == 1:
        # Open the file for writing
        with open('results.txt', 'w') as file:
            # Writing the initial condition
            file.write('PARAMETERS obtained by calibration:\n')
            file.write('X0 used:\n')
            file.write('\n'.join(map(str, initial_cond)))
            file.write('\n-----------------------\n')

            # Writing the USA market parameters
            file.write('Calibrated parameters for the USA market:\n')
            file.write('kappa_USA: {}\n'.format(kappa_USA))
            file.write('theta_USA: {}\n'.format(theta_USA))
            file.write('sigma_USA: {}\n'.format(sigma_USA))
            file.write('-----------------------\n')

            # Writing the EU market parameters
            file.write('Calibrated parameters for the EU market:\n')
            file.write('kappa_EU: {}\n'.format(kappa_EU))
            file.write('theta_EU: {}\n'.format(theta_EU))
            file.write('sigma_EU: {}\n'.format(sigma_EU))
            file.write('-----------------------\n')

    # Display initial condition
    print('-----------------------')
    print('X0 used:')
    for val in initial_cond:
        print(val)
    print('PARAMETERS obtained by calibration:')
    print('-----------------------')
    
    # Display calibrated parameters for the USA market
    print('Calibrated parameters for the USA market:')
    print('kappa_USA:', kappa_USA)
    print('theta_USA:', theta_USA)
    print('sigma_USA:', sigma_USA)
    print('-----------------------')
    
    # Display calibrated parameters for the EU market
    print('Calibrated parameters for the EU market:')
    print('kappa_EU:', kappa_EU)
    print('theta_EU:', theta_EU)
    print('sigma_EU:', sigma_EU)
    print('-----------------------')

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

def blk_calibration(sigma, dataset, F0, B0, date_settlement):
    # Conventions
    conv_ACT365 = 3

    # Initial parameters
    MSE = np.zeros(len(dataset['datesExpiry']))

    # Computation
    for ii in range(len(dataset['datesExpiry'])):
        # Initialization
        put_length = len(dataset['putAsk'][ii]['prices'])

        # Time to maturity
        TTM = yearfrac(date_settlement, datenum(dataset['datesExpiry'][ii]), conv_ACT365)
        
        # Interest rate
        interest_rate = -np.log(B0[ii]) / TTM

        # Pricing
        # Price of Call/Puts through the Black formula
        prices = blkprice(F0[ii], dataset['strikes'][ii]['value'], interest_rate, TTM, sigma)
        
        call_prices = prices[put_length:]
        put_prices = prices[:put_length] - F0[ii] * B0[ii] + dataset['strikes'][ii]['value'][:put_length] * B0[ii]
        
        mean_call_price = (dataset['callAsk'][ii]['prices'] + dataset['callBid'][ii]['prices']) / 2
        mean_put_price = (dataset['putAsk'][ii]['prices'] + dataset['putBid'][ii]['prices']) / 2

        # Computation of RMSE
        MSE[ii] = np.sqrt(np.mean((np.concatenate([call_prices, put_prices]) - np.concatenate([mean_call_price, mean_put_price])) ** 2))

    # Final adjusting of RMSE
    return np.sum(MSE)


# Function to interpolate interest rates
def rate_interpolation(dates, B0, settlement_date, interp_date):
    TTM = np.array([(date - settlement_date).days / 365 for date in dates])
    rate = -np.log(B0) / TTM
    rate_interpolator = interp1d(TTM, rate)
    return rate_interpolator((interp_date - settlement_date).days / 365)

# Function to extract interest rate and time to maturity
def interp_pricing_params(dates, B0, date_settlement, year_to_maturity):
    conv_ACT365 = 365  # Assuming ACT/365 convention

    # Calculate interpolation date
    interp_date = (datetime.strptime(date_settlement, '%Y-%m-%d') - timedelta(days=1) + timedelta(days=365 * year_to_maturity)).date()

    # Calculate time to maturity
    TTM = (interp_date - datetime.strptime(date_settlement, '%Y-%m-%d').date()).days / conv_ACT365

    # Calculate interest rate
    rate = rate_interpolation(dates, B0, datetime.strptime(date_settlement, '%Y-%m-%d').date(), interp_date)

    return rate, TTM

import numpy as np

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
    nSim = 1000000
    
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
    
    Y_1 = -Beta_USA * gamma_USA**2 * G_1 * TTM + gamma_USA * np.sqrt(TTM * G_1) * g_1
    Y_2 = -Beta_EU * gamma_EU**2 * G_2 * TTM + gamma_EU * np.sqrt(TTM * G_2) * g_2
    Z = -Beta_z * gamma_z**2 * G_z * TTM + gamma_z * np.sqrt(TTM * G_z) * g_z
    
    # Marginal processes
    X_1 = Y_1 + a_USA * Z
    X_2 = Y_2 + a_EU * Z
    
    Xt = np.column_stack((X_1, X_2))
    
    stock = S0 * np.exp((rates + drift_compensator) * TTM + Xt)
    
    return stock, S0

import numpy as np

def stock_simulation_Black(sigmas, F0, rates, rho, TTM):
    nSim = int(1e6)

    # Simulation of the NIG process
    covarianceMatrix = np.array([[TTM, rho * TTM], [rho * TTM, TTM]])
    cholesky_matrix = np.linalg.cholesky(covarianceMatrix)

    meanVector = np.array([0, 0])

    np.random.seed(2)
    g = np.random.multivariate_normal(meanVector, covarianceMatrix, nSim)

    # Creation of Xt dynamic
    Xt = -0.5 * sigmas ** 2 * TTM + sigmas * np.sqrt(TTM) * g.T

    # Computation of the initial stock
    prices = F0 * np.exp(Xt)

    # Computation for the antithetic behavior
    Xt_AV = -0.5 * sigmas ** 2 * TTM - sigmas * np.sqrt(TTM) * g.T
    pricesAV = F0 * np.exp(Xt_AV)

    return prices, pricesAV

def blk_semiclosed(s1_0, rate1, rate2, sigma1, sigma2, rho, TTM):
    discount = np.exp(-rate1 * TTM)

    def d1(w):
        return ((rate1 - sigma1 ** 2 / 2) * TTM + rho * w * sigma1 + (1 - rho ** 2) * sigma1 ** 2 * TTM) / (sigma1 * np.sqrt((1 - rho ** 2) * TTM))

    def d2(w):
        return d1(w) - np.sqrt((1 - rho ** 2) * TTM) * sigma1

    def integrand(w):
        return (np.exp(rate1 * TTM - sigma1 ** 2 * rho ** 2 * TTM / 2 + sigma1 * rho * w) * norm.cdf(d1(w)) - norm.cdf(d2(w))) * np.exp(-w ** 2 / (2 * TTM)) / np.sqrt(2 * np.pi * TTM)

    xmax = (np.log(0.95) - (rate2 - sigma2 ** 2 / 2) * TTM) / sigma2

    price, _ = quad(integrand, -np.inf, xmax)
    price *= s1_0 * discount

    return price

def disp_contract_prices(price_Levy, CI_Levy, price_Blk, CI_Blk, price_Blk_AV, CI_Blk_AV, price_SemiclosedBlk):
    print("Prices of the derivative")
    print("-----------------------------------------------------------------------------------------------------")
    print("| Model              |  Price               |  Confidence Interval                       |")
    print("-----------------------------------------------------------------------------------------------------")
    print("| Levy:              |  {:.8f}%    |  [ {:.8f}% ,  {:.8f}% ]".format(price_Levy, CI_Levy[0], CI_Levy[1]))
    print("| Black:             |  {:.8f}%    |  [ {:.8f}% ,  {:.8f}% ]".format(price_Blk, CI_Blk[0], CI_Blk[1]))
    print("| Black AV:          |  {:.8f}%    |  [ {:.8f}% ,  {:.8f}% ]".format(price_Blk_AV, CI_Blk_AV[0], CI_Blk_AV[1]))
    print("| Black semi-closed: |  {:.8f}%    |  --".format(price_SemiclosedBlk))
    print("-----------------------------------------------------------------------------------------------------")
