clear

data = xlsread('usdata_ypr.xls');

T = size(data,1);
ncol = size(data,2);

data = data(2:end,2:ncol); % first column is time 

chi = [0.0045;1.005;1.005;0.75;1.7;0.3;0.5;0.5;0.5;0.5;0.5;0.5;0.5];

% lhout = likelihood(chi,data);
% 
% logpostout = logpostwout(chi,data);
% 
% maxpostout = maxpost(chi,data);

finalqn(chi,date);