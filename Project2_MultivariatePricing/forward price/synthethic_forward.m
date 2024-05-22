function [Gi, Gi_ask, Gi_bid] = synthethic_forward(call_bid, call_ask, put_bid, put_ask)
% Computation of the Synthetic Forward Gi for each of the required
% horizons, the SF is the building block of the Forward computation
% 
% INPUT:
% call_bid:          [VEC] bid prices of the Call
% call_ask:          [VEC] ask prices of the Call
% put_bid:           [VEC] bid prices of the Put
% put_ask:           [VEC] ask prices of the Put
% 
% OUTPUT:
% Gi:                Synthetic Forward for each strike
% Gi_ask:            Synthetic Forward Ask for each strike
% Gi_bid:            Synthetic Forward Bid for each strike

    %% Computation Synthetic Forward Bid

    Gi_bid = call_bid - put_ask;

    %% Computation Synthetic Forward Ask

    Gi_ask = call_ask - put_bid;

    %% Computation final Synthetic Forward

    Gi = (Gi_bid + Gi_ask)/2;

end % function Synthethic_Forward