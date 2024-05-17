function [Gi, Gi_ask, Gi_bid] = Synthethic_Forward(call_bid, call_ask, put_bid, put_ask, index)
% Computation of the Synthetic Forward Gi for each of the required
% horizons, the SF is the building block of the Forward computation
% 
% INPUT:
% call_bid:          [STRUCT] bid prices of the Call
% call_ask:          [STRUCT] ask prices of the Call
% put_bid:           [STRUCT] bid prices of the Put
% put_ask:           [STRUCT] ask prices of the Put
% index:             [SCALAR] index related to the maturity considered
% 
% OUTPUT:
% Gi:                Synthetic Forward for each strike

    %% Computation Synthetic Forward Bid
    
    % Computation of the necessary prices from CallBid and PutAsk
    call_bid_prices = call_bid(index).prices;
    put_ask_prices = put_ask(index).prices;

    Gi_bid = call_bid_prices - put_ask_prices;

    %% Computation Synthetic Forward Ask

    % Computation of the necessary prices from CallAsk and PutBid
    call_ask_prices = call_ask(index).prices;
    put_bid_prices = put_bid(index).prices;

    Gi_ask = call_ask_prices - put_bid_prices;

    %% Computation final Synthetic Forward

    Gi = (Gi_bid + Gi_ask)/2;

end % function Synthethic_Forward