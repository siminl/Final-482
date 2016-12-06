% ~~~~ Final Exam Questions ~~~~ %

function rt = finalqn(data)
   chi = [0.0045;1.005;1.005;0.75;1.7;0.3;0.5;0.5;0.5;0.5;0.5;0.5;0.5];
 
   postout = maxpost(chi,data);
   n = size(data,2); % should be 3
   
   % initialize MCMC with the posterior modes
   chi = postout.chimode; % chi mode
   sighat = postout.ih; % inv hessian mode
   display(sighat)
   post = postout.postmode;
   lh = postout.lhmode;
   prior = postout.priormode;
   H = postout.Hmode;
   M = postout.Mmode;
   G1 = postout.G1mode;
   
   %% Metropolis MCMC
    iter = 5000;
    horizon = 30;
    postchi = zeros(iter,length(chi));
    
    irf = zeros(horizon,n,iter);
    arcount = 0; % acceptance rate tracker
    trplots = zeros(iter,length(chi));
    for i = 1:iter
        k = 2;
        V = (k^2)* sighat; 

        % new draw
        chitilde = mvnrnd(chi,V);    

%         chitilde = mvnrnd(zeros(length(chi),1),V);            
        logposttilde = logpostwout(chitilde,data);   
        disp(logposttilde);
        if logposttilde.lh ~= 0
            posttilde = -logposttilde.postprob;
            Htilde = logposttilde.H;
            Mtilde = logposttilde.M;
            G1tilde = logposttilde.G1;
        else
            posttilde = 0;
        end

        % discard the ones with eu showing nonuniqueness or nonexistence
        % note that in the likeli hood function, I have set likelihood = 0 if eu
        % ~=[1 1]; therefore, alpha = 0 below, the new draw is rejected

        alpha = min(1,posttilde/post);
        disp(alpha)

        % Acceptance-rejection step
        pAR = unifrnd(0,1,1,1);
        if pAR <= alpha
            chi = chitilde;
            H = Htilde;
            M = Mtilde;
            G1 = G1tilde;
            post = posttilde;
            
            arcount = arcount +1;
        else 
            chi = chi;
            H = H;
            M = M;
            G1 = G1;
            post = post; 
            
        end;
        
        
        postchi(i,:) = chi; % store chi for part (e)&(g)
        %% ~~~~~~~~~~~~~~~~ irf ~~~~~~~~~~~~~~~~~~~ %%

        % From calculation: IRF = H*G^s(chi)*M(chi)*A
        % where A = diag(sigMP,sigz,sigtheta,sigb) transform var(epsilon) into an
        % identity matrix

        % form the A matrix 
        A = diag(chi(end-3:end));

        impulse = zeros(size(A,2),1); 
        impulse(1,1) = 1; % the impulse in monetary shock

        for h = 1:horizon         
            temp = H* G1^h * M*A* impulse;
            irf(h,:,i) = temp(1:n,:);
        end;


        % trace plots to investigate simulation
        trplots(i,:) = chi';

    end; % end of MH MCMC
    
%% ~~~~~~ trace plots ~~~~~~~ %%
    figure 
    subplot(2,3,1)
    plot(trplots(:,1));
    subplot(2,3,2)
    plot(trplots(:,3));
    subplot(2,3,3)
    plot(trplots(:,5));
    subplot(2,3,4)
    plot(trplots(:,7));
    subplot(2,3,5)
    plot(trplots(:,10));
    subplot(2,3,6)
    plot(trplots(:,14));

    annotation('textbox', [0 0.9 1 0.1], ...
        'String', 'Trace Plots', ...
        'EdgeColor', 'none', ...
        'FontSize',14,...
        'Color','black',...
        'HorizontalAlignment', 'center')

%% ~~~~~~~~~~~~~~~~~ part (f) IRF error bands ~~~~~~~~~~~~~~~ %%

    eb = zeros(horizon,n,4);

    % calculate IRF at the modes
    chimode = postout.xh; % chi mode
    Hmode = postout.Hmode;
    Mmode = postout.Mmode;
    G1mode = postout.G1mode;


    % form the A matrix 

    Amode = diag(chimode(end-3:end));

    impulse = zeros(size(Amode,2),1); 
    impulse(1,1) = 1; % the impulse in monetary policy shock

    for h = 1:horizon         
        temp = Hmode* G1mode^h * Mmode*Amode* impulse;
        eb(h,:,4) = temp(1:n,:);
    end;

    trim = 0.2;
    
    % calculate median and 20 quarters

    for i=1:horizon
        for j = 1:n
            eb(i,j,1:3) = quantile(irf(i,j,iter*trim:iter),[0.1 0.9 0.5]);   
        end;
    end;

    rt.eb = eb;

    %~~ plot error bands ~~%
    hrz = 1:horizon;
    
    figure 
    subplot(3,1,1)  
    plot(hrz, eb(:,1,1),'--m',hrz,eb(:,1,2),'--m',hrz,eb(:,1,3),'-k',hrz,eb(:,1,4),'-b')
    title('Log-Output')
 
    subplot(3,1,2) 
    plot(hrz, eb(:,2,1),'--m',hrz,eb(:,2,2),'--m',hrz,eb(:,2,3),'-k',hrz,eb(:,1,4),'-b')
    title('Inflation Rate')
    legend('Lower Bound','Upper Bound','Median', 'Mode')

    subplot(3,1,3) 
    plot(hrz, eb(:,3,1),'--m',hrz,eb(:,3,2),'--m',hrz,eb(:,3,3),'-k',hrz,eb(:,1,4),'-b')
    title('Interest Rate')

    annotation('textbox', [0 0.9 1 0.1], ...
        'String', 'IRF with 90% Error Bands', ...
        'EdgeColor', 'none', ...
        'FontSize',14,...
        'Color','black',...
        'HorizontalAlignment', 'center')

    
 %% ~~~~~~~~~ part (g) ML ~~~~~~~~~~~ %%
 
 %% ~~~~~~~~~ g.1 second order approximation ~~~~~~~~~ %%
   postmode = postout.postmode;
   logsighatmode = postout.logsighat;
   
   MLsecorder = postmode* T^(-n)* det(logsighatmode)* (2*pi)^(n/2);
   
 %% ~~~~~~~~~ g.2 Geweke's Modified Harmonic Mean ~~~~~~~~~ %%
 
 chihat = mean(postchi,1); % mean over iteration, check
 
 lam = [0.3;0.5;0.9];
 MLMH = zeros(length(lam),1);
 
 for i = 1:length(lam)
     fi = rowfun(@GewekeMH, postchi,'chihat',chihat,'logsighatmode',logsighatmode,'lambda',lam(i),'n',n);
     fbar = mean(fi,1); % mean over iterations, check
     logMLMHinv = fbar + mean(fi-fbar,1);
     MLMH(i) = exp(-logMLMHinv);
 end
 

rt.MLsecorder = MLsecorder;
rt.MLMH = MLMH;
rt.AR = arcount/iter; % report acceptance rate
end

% ~~~~ MCMC ~~~~ %
