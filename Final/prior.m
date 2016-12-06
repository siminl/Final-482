% ~~~~ prior distribution ~~~~ %

function logpriorprob = prior(chi)
% chi = (gamma,pi,r,xip,phipi,phiy,rhoR,rhoz,rhob,...
%    sigMP,sigz,sigtheta,sigb) is the parameter vector
% nu is a constant 

gamma = chi(1);
pi = chi(2);
r = chi(3);

xi = chi(4);
phipi = chi(5);
phiy = chi(6);
rhoR = chi(7);
rhoz = chi(8);
rhob = chi(9);
sigMP = chi(10);
sigz = chi(11);
sigtheta = chi(12);
sigb = chi(13);

% normal
gammamu = 0.0045;
gammasd = 0.00025;
pimu = 1.005;
pisd = 0.0015;
rmu = 1.005;
rsd = 0.0015;
phipimu = 1.7;
phipisd = 0.3;
phiymu = 0.3;
phiysd = 0.15;

nu = 2; % constant

%beta
xia = 5.5; 
xib = 1.8333; 
rhoRa = 1.5;
rhoRb = 1.5;
rhoza = 1.5;
rhozb = 1.5;
rhoba = 1.5;
rhobb = 1.5;

% IG-1
sigMPd = 2.0395;
sigMPs = 0.1679;
sigzd = 2.0395;
sigzs = 0.1679;
sigthetad = 2.0395;
sigthetas = 0.1679;
sigbd = 2.0395;
sigbs = 0.1679;


logpriorprob = log(normpdf(gamma,gammamu,gammasd))+log(normpdf(pi,pimu,pisd))+log(normpdf(r,rmu,rsd))+...
    log(nu)+log(betapdf(xi,xia,xib))+log(normpdf(phipi,phipimu,phipisd))+log(normpdf(phiy,phiymu,phiysd))+...
    log(betapdf(rhoR,rhoRa,rhoRb))+log(betapdf(rhoz,rhoza,rhozb))+log(betapdf(rhob,rhoba,rhobb))+...
    log(ig1pdf(sigMP,sigMPd,sigMPs))+log(ig1pdf(sigz,sigzd,sigzs))+log(ig1pdf(sigtheta,sigthetad,sigthetas))+...
    log(ig1pdf(sigb,sigbd,sigbs));


end