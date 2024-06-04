function price = blk_semiclosed(s1_0, rate1, rate2, sigma1, sigma2, rho, TTM)
%
% INPUT:
% s1_0:         [SCALAR] initial stock value
% rate1         [SCALAR] rate USA mkt
% rate2         [SCALAR] rate EU mkt
% sigma1        [SCALAR] sigma USA mkt
% sigma2        [SCALAR] sigma EU mkt
% rho           [SCALAR] historical correlation
% TTM           [SCALAR] time to maturity
% OUTPUT:
% price : price od the option obtained via the semiclosed formula
% USES : blk_semiclosed()


    % discount = exp(-rate*TTM);
    % 
    % d1 = @(w) ((rate-sigma1^2/2)/sigma1 + (1-rho^2)*sigma1 + rho * w) / sqrt(1-rho^2);
    % d2 = @(w) d1(w) - sqrt(1-rho^2)*sigma1;
    % 
    % integrand = @(w) ( exp(rate-sigma1^2*rho^2/2+sigma1*rho*w) .* normcdf(d1(w)) - normcdf(d2(w)) ) .* exp(-w.^2./2)./sqrt(2*pi);
    % 
    % xmax = (log(0.95)-(rate-sigma2^2/2))/sigma2;
    % 
    % price = s1_0*discount* integral( integrand, -inf, xmax );


    discount = exp(-rate1*TTM);
    
    d1 = @(w) ((rate1-sigma1^2/2)*TTM + rho*w*sigma1 + (1-rho^2)*sigma1^2*TTM ) / (sigma1*sqrt((1-rho^2)*TTM));
    d2 = @(w) d1(w) - sqrt((1-rho^2)*TTM)*sigma1;
    
    integrand = @(w) ( exp(rate1*TTM-sigma1^2*rho^2*TTM/2+sigma1*rho*w) .* normcdf(d1(w)) - normcdf(d2(w)) ) .* exp(-w.^2./(2*TTM))./sqrt(2*pi*TTM);
    
    xmax = (log(0.95)-(rate2-sigma2^2/2)*TTM)/sigma2;
    
    price = s1_0*discount* integral( integrand, -inf, xmax );
end