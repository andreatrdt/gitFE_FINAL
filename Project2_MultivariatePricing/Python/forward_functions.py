# Built in libraries

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.io as sio
import datetime
import warnings
import os

'''
from scipy.stats import norm
from dateutil.relativedelta import relativedelta
from scipy.integrate import quad
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
from scipy.optimize import OptimizeResult
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from scipy.optimize import brentq '''


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
    
    # Mean values of the strikes and synthetic forwards
    G_hat = np.mean(Gi)
    K_hat = np.mean(Ki)
    
    # Ensure the shapes are correct for dot product
    Ki_diff = Ki - K_hat
    Gi_diff = Gi - G_hat
    
    # Final computation of B_bar(0, T)
    B_bar = -np.dot(Ki_diff.flatten(), Gi_diff.flatten()) / np.dot(Ki_diff.flatten(), Ki_diff.flatten())
    
    return B_bar


# Syntehtic forward
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
def forward_prices(dataset):
    """
    Computation of the forward prices following the Baviera, Azzone paper
    
    INPUT:
    dataset:           [CLASS] data containing all the required tables
    flag:              [0: without plots & slope; 1: with plots & slope]
    
    OUTPUT:
    F_vector:          [ARRAY] forwards value F(0, T)
    B_bar_vector:      [ARRAY] market calibrated discount B(0, T)
    
    USES:
    function synthetic_forward()
    function estimation_discount_factor()
    """

    # Inizialization
    length_dataset = dataset._length_dates()

    F_vector = np.zeros(length_dataset)
    B_bar_vector = np.zeros(length_dataset)

    # Computation of the forwards
    for ii in range(length_dataset):

        # Computation of the strikes and synthetic forwards
        Ki = dataset._strikes_i(ii)

        # Synthetic forwards
        Gi, Gi_ask, Gi_bid = synthetic_forward(dataset._callBid_i(ii), dataset._callAsk_i(ii),
                                               dataset._putBid_i(ii), dataset._putAsk_i(ii))
        
        # Computation of the estimated discount factor
        B_bar_vector[ii] = estimation_discount_factor(Gi, Ki)

        # Computation of the forward prices
        F_i_vector = Gi / B_bar_vector[ii] + Ki
        F_ask__i_vector = Gi_ask / B_bar_vector[ii] + Ki
        F_bid__i_vector = Gi_bid / B_bar_vector[ii] + Ki

        # Computation of the required Forwards
        F_vector[ii] = np.mean(F_i_vector)
        

    return F_vector, B_bar_vector