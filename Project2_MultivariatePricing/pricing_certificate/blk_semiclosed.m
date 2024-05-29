function price = blk_semiclosed(s1_0, rate, sigma1, sigma2, rho, TTM)

discount = exp(-rate*TTM);

d1 = @(w) ((rate-sigma1^2/2)/sigma1 + (1-rho^2)*sigma1 + rho * w) / sqrt(1-rho^2);
d2 = @(w) d1(w) - sqrt(1-rho^2)*sigma1;

integrand = @(w) ( exp(rate-sigma1^2*rho^2/2+sigma1*rho*w) .* normcdf(d1(w)) - normcdf(d2(w)) ) .* exp(-w.^2./2)./sqrt(2*pi);

xmax = (log(0.95)-(rate-sigma2^2/2))/sigma2;

price = s1_0*discount* integral( integrand, -inf, xmax );



end