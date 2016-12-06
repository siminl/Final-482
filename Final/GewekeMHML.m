%% ~~~~ intermediate step in calculating Geweke's ML ~~~~ %%
%% ~~~~ calculate f_i(chi) ~~~~ %%
function rt = GewekeMHML(chirow,chihat,logsighatmode,lambda,n)
    secorder = (chirow-chihat)'*inv(logsighatmode)*(chirow-chihat);
    rt = (2*pi)^(-n/2)*det(logsighatmode)^(-1/2)*exp(-1/2*secorder)* ...
        (secorder <= chi2inv(1-lambda,n));
end