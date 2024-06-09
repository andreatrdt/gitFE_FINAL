# Data set class for the option data
#
# Authors :
# M. Maspes
# A. Tarditi
# M. Torba

#%%

# Import of the required libraries
import numpy as np
import pandas as pd


class Data_maturity:
    def __init__(self, data_structure, idx):

        # Initialization of the structure to a class object
        self._callAsk = np.squeeze(data_structure[0, 0]['callAsk'][0, idx]['prices'])
        self._callBid = np.squeeze(data_structure[0, 0]['callBid'][0, idx]['prices'])

        self._putAsk = np.squeeze(data_structure[0, 0]['putAsk'][0, idx]['prices'])
        self._putBid = np.squeeze(data_structure[0, 0]['putBid'][0, idx]['prices'])

        self._strikes = np.squeeze(data_structure[0, 0]['strikes'][0, idx]['value'])

    # Getters:
    def _callAsk(self):
        return self._callAsk

    def _callBid(self):
        return self._callBid

    def _putAsk(self):
        return self._putAsk

    def _putBid(self):
        return self._putBid

    def _strikes(self):
        return self._strikes

    # Setters
    def _callAsk(self, new_value):
        self._callAsk = new_value

    def _callBid(self, new_value):
        self._callBid = new_value

    def _putAsk(self, new_value):
        self._putAsk = new_value

    def _putBid(self, new_value):
        self._putBid = new_value

    def _strikes(self, new_value):
        self._strikes = new_value


class Dataset:

    def __init__(self, data_structure):
        # Initialization of the structure to a class object

        # Initialize an empty vector
        self._list_obj = np.zeros(np.shape(data_structure[0, 0]['callAsk'])[1], dtype=Data_maturity)
        # self._list_dates = np.zeros(np.shape(data_structure[0, 0]['callAsk'])[1], dtype=DateTime)
        self._list_dates = []

        self._spot = np.squeeze(data_structure[0, 0]['spot'][0, 0])

        # Populate the list
        for ii in range(np.shape(data_structure[0, 0]['callAsk'])[1]):

            # Create object for the ith maturity
            obj_data = Data_maturity(data_structure, int(ii))

            self._list_obj[ii] = obj_data


    # Setters
    def _list_dates(self, new_value):
        self._list_dates = new_value

    def _set_callAsk(self, idx, new_value):
        self._list_obj[idx]._callAsk = new_value

    def _set_callBid(self, idx, new_value):
        self._list_obj[idx]._callBid = new_value

    def _set_putAsk(self, idx, new_value):
        self._list_obj[idx]._putAsk = new_value

    def _set_putBid(self, idx, new_value):
        self._list_obj[idx]._putBid = new_value

    def _set_strikes(self, idx, new_value):
        self._list_obj[idx]._strikes = new_value

    # Getters
    def _spot(self):
        return self._spot

    def _list_dates(self):
        return self._list_dates

    def _list_dates_i(self, idx):
        return self._list_dates[idx]

    def _length_dates(self):
        return len(self._list_dates)

    def _strikes_i(self, idx):
        return self._list_obj[idx]._strikes

    def _callBid_i(self, idx):
        return self._list_obj[idx]._callBid

    def _callAsk_i(self, idx):
        return self._list_obj[idx]._callAsk

    def _putBid_i(self, idx):
        return self._list_obj[idx]._putBid

    def _putAsk_i(self, idx):
        return self._list_obj[idx]._putAsk

    def _remove_last_element(self):
        self._list_obj = self._list_obj[:-1]