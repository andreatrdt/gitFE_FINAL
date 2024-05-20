function rho = hist_corr(data)

returns_EU = data.Annually(:,1);

returns_USA = data.Annually(:,2);

rho = corr(returns_EU,returns_USA);
end