function prob = ig1pdf(sig,d,s)

C = (gamma(d/2)*(2/s)^(d/2))/2;
prob = C^(-1) *sig^(-(d+1))*exp(-s/(2*(sig^2)));

end