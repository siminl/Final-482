% ~~~~ likelihood function ~~~~ %

function lh = likelihood(chi, data)

 T = size(data,1);
 n = size(data,2);
 zetalen = 8;
    % parameter vector
    gamma = chi(1);
    pi = chi(2);
    r = chi(3);
    nu = 2;
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

    % data: the observed dataset, which will be used in the kalman filter
   
    % the steady state solutions 
    beta = exp(gamma)/r; % the discount faction 
    R = pi*r;

    % Form Gensys canonical form matrix
    g00 = [1 -(1-beta*xi)*(1-xi)*(1+nu)/((1+beta)*xi) 0 0 0 -beta/(1+beta) 0 0;
        0 1 1 -rhoz -(1-rhob) -1 -1 0; 
        -(1-rhoR)*phipi -(1-rhoR)*phiy 1 -(1-rhoR)*phiy 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 1 0 0 0;
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1];
    
    g11 = [1/(1+beta) 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 -(1-rhoR)*phiy rhoR 0 0 0 0 0;
        0 0 0 rhoz 0 0 0 0;
        0 0 0 0 rhob 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        0 1 0 0 0 0 0 0];
    
    psi = [0 0 1 0;0 0 0 0;1 0 0 0;0 1 0 0;...
        0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0];
    
    c = zeros(8,1); % length of zeta = 8
    
    Pigensys = zeros(zetalen,2); % length of zeta = 8
    Pigensys(6,1) = 1;
    Pigensys(7,2) = 1;
    
    % call Gensys
    [G1,C,impact,fmat,fwt,ywt,gev,eu,loose] = gensys(g00,g11,c,psi,Pigensys);
    % output model 
    % y(t)=G1*y(t-1)+impact*z(t)+C
    
    % sanity check; discard the output if output chi not unique or not
    % exist.
        if eu(1) == 1 || eu(1) == 1
                
        % the observational function
        H = [0 1 0 1 0 0 0 -1;4 0 0 0 0 0 0 0;0 0 4 0 0 0 0 0];
        Cchi = [100*gamma;400*(pi-1);400*(R-1)];

        % Kalman filter

        % initialize the KF
        shat = zeros(zetalen,1); % length of zeta_t = 8;
        % var(epsilon_t) = I, since the fundamental shocks are indept
        Q = impact*diag([sigMP^2,sigz^2,sigtheta^2,sigb^2])*impact'; % V(z(t)) = Q = impact*diag(sigMP,sigz,sigtheta,sigb)*impact'
        sig = dlyap(G1,Q); % solution to sig = G1(chi)*sig*G1(chi)'+impact(chi)*var(epsilon_t)*impact(chi)'

        Robserr = zeros(n,n);        
        omega=G1*sig*G1'+Q;
        sigma=H*omega*H'+Robserr;

        % P(y_1|chi)
        loglh=-.5*log(det(sigma))-.5*((data(1,:)'-Cchi-H*G1*shat)'/sigma)*(data(1,:)'-Cchi-H*G1*shat)-.5*n*log(2*pi);
        
        sigmat = zeros(zetalen,zetalen,T); % store posterior variance        
        sigmat (:,:,1) = sig;
        
        % iterations 

        for i = 2:T
            [shatnew,signew,loglhtemp]=kfilter(data(i,:),H,Cchi,G1,shat,sig,Robserr,Q);
            shat = shatnew;
            sig = signew;
            loglh = loglh+loglhtemp;
            sigmat(:,:,i) = sig;
        end
        
        lh.sigmat = sigmat;
        
        lh.loglh = loglh;        
        lh.G1 = G1;
        lh.M = impact;
        lh.Q = Q;
        lh.H = H; % observational function
        
    else
        % assign dummies
        lh.loglh = 0;
        lh.G1 = [];
        lh.M = [];
        lh.Q = [];
        lh.H = []; % observational function
        lh.sigmat = [];
    end;

end