function [Gi, Gi_ask, Gi_bid] = synthethic_forward(call_bid, call_ask, put_bid, put_ask)
% Computation of the Synthetic Forward Gi for each of the required
% horizons, the SF is the building block of the Forward computation
% 
% INPUT:
% call_bid:          [VECTOR] bid prices of the Call
% call_ask:          [VECTOR] ask prices of the Call
% put_bid:           [VECTOR] bid prices of the Put
% put_ask:           [VECTOR] ask prices of the Put
% 
% OUTPUT:
% Gi:                [VECTOR] Synthetic Forward for each strike
% Gi_ask:            [VECTOR] Synthetic Forward Ask for each strike
% Gi_bid:            [VECTOR] Synthetic Forward Bid for each strike
%
% USES: none

% Authors:
% M.Maspes, A.Tarditi, M.Torba

    %% Computation Synthetic Forward Bid

    Gi_bid = call_bid - put_ask;

    %% Computation Synthetic Forward Ask

    Gi_ask = call_ask - put_bid;

    %% Computation final Synthetic Forward

    Gi = (Gi_bid + Gi_ask)/2;

end % function Synthethic_Forward