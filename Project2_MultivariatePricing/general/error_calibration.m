function [error_call_prices, error_put_prices]=error_calibration(model_call_prices, model_put_prices, mkt_call_prices_bid, mkt_call_prices_ask, mkt_put_prices_bid, mkt_put_prices_ask)
% Computation of the model error in the pricing of put and call options
% with respect to market bid and ask prices for a fixed expiry.
% 
% INPUT:
% model_call_prices:    [VECTOR] call prices obtained via the model for the
%                       considered expiry
% model_put_prices:     [VECTOR] put prices obtained via the model for the
%                       considered expiry
% mkt_call_prices_bid:  [VECTOR] call prices of the bid taken from the
%                       market for the considered expiry
% mkt_call_prices_ask:  [VECTOR] call prices of the ask taken from the
%                       market for the considered expiry
% mkt_put_prices_bid:   [VECTOR] put prices of the bid taken from the
%                       market for the considered expiry
% mkt_put_prices_ask:   [VECTOR] put prices of the ask taken from the
%                       market for the considered expiry
% 
% OUTPUT:
% error_call_prices:    [VECTOR] error in the pricing of the call options
% error_put_prices:     [VECTOR] error in the pricing of the put options
%
% USES:     none

% Authors:
% M.Maspes, A.Tarditi, M.Torba


%% Call prices error

% Indices to find when the call model price is not in the bid-ask interval:
idx1_call =  model_call_prices > mkt_call_prices_ask;
idx2_call = model_call_prices < mkt_call_prices_bid;

% Computation of the error on the call prices:
error_call_prices = 100*((model_call_prices-mkt_call_prices_ask)./mkt_call_prices_ask.*idx1_call + (mkt_call_prices_bid-model_call_prices)./mkt_call_prices_bid.*idx2_call);

%% Put prices error

% Indices to find when the put model price is not in the bid-ask interval:
idx1_put =  model_put_prices > mkt_put_prices_ask;
idx2_put = model_put_prices < mkt_put_prices_bid;

% Computation of the error on the put prices:
error_put_prices = 100*((model_put_prices-mkt_put_prices_ask)./mkt_put_prices_ask.*idx1_put + (mkt_put_prices_bid-model_put_prices)./mkt_put_prices_bid.*idx2_put);

end