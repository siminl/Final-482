% calculate the log posterior 

function postprob = logpost(chi, data)

    logpriorprob = prior(chi);

    lh = likelihood(chi, data);
    postprob = -(logpriorprob+ lh.loglh); 

%     rt.postprob = postprob;
%     
%     rt.lh = lh.loglh;
%     rt.G1 = lh.G1;
%     rt.M = lh.M;
%     rt.H = lh.H;
%     rt.sigmat = lh.sigmat;
%     rt.logsighat = lh.logsighat;
%     
%     rt.logpriorprob = logpriorprob;
    
end