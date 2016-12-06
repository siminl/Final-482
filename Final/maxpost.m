% ~~~~ max log posterior using maximization routine ~~~~ %

function rt = maxpost(chi,data)

%  Hini = diag(chi);

 Hini = 0.1*eye(13,13);
% minimization routine to get argmax chi
[fh,xh,gh,H,itct,fcount,retcodeh] = ...
    csminwel(@logpost, chi, Hini,[],1e-6,10^20,data);

rt.chimode = xh; % argmax chi
rt.ih = H; % inverse Hessian 

logpostout = logpostwout(xh,data);

% postprob = logpostout.postprob;
% lh = logpostout.lh;
% priorprob = logpostout.priorprob;

rt.postmode = -logpostout.postprob;
rt.lhmode = logpostout.lh;
rt.priormode = logpostout.logpriorprob;
rt.Hmode = logpostout.H;
rt.Mmode = logpostout.M;
rt.G1mode = logpostout.G1;
rt.Qmode = logpostout.Q;


end
