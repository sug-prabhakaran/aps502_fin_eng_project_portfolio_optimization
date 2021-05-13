%% Appendix - MATLAB Code & Results Output
%---------------------------------------------------------------
% PROBLEM 1 - PART 1 
%---------------------------------------------------------------
% As per matrix formulation in report, create vectors and matrices:

% obj function coefficients
f = [108; 94; 99; 92.7; 96.6; 95.9; 92.9; 110; 
     104; 101; 107; 102; 95.2; 0; 0; 0; 0; 0];

% inequality constraint matrix A
A = [-10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -10 -7 -100 1 0 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -110 -107 0 -1.02 1 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -110 -108 -106 0 0 0 0 -1.03 1 0 0;
     -10 -7 -8 -6 -7 -106 -105 0 0 0 0 0 0 0 0 -1.04 1 0;
     -10 -7 -8 -106 -107 0 0 0 0 0 0 0 0 0 0 0 -1.05 1; 
     -100 -107 -108 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.06];

% inequality constraint vector b
b = [-500; -200; -800; -400; -700; -900];

% lower bound and upper bound vectors
lb = zeros(18,1);
ub = inf*ones(18,1);

% calculate objective function value fval and decision var x
[x, fval] = linprog(f,A,b,[],[],lb,ub)

%--------------------------------------------------------------
% PROBLEM 1 - PART 2
%--------------------------------------------------------------
% updated constraint matrix A2 with <=50% B-rated bond row
A2 = [-10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -10 -7 -100 1 0 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -110 -107 0 -1.02 1 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -110 -108 -106 0 0 0 0 -1.03 1 0 0;
     -10 -7 -8 -6 -7 -106 -105 0 0 0 0 0 0 0 0 -1.04 1 0;
     -10 -7 -8 -106 -107 0 0 0 0 0 0 0 0 0 0 0 -1.05 1; 
     -100 -107 -108 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.06;
     108 94 99 92.7 96.6 95.9 -92.9 -110 -104 
     -101 -107 -102 -95.2 0 0 0 0 0];

% updated constraint vector b2
b2 = [-500; -200; -800; -400; -700; -900; 0];

% calculate new objective function value fval2 and dec. var x2
[x2, fval2] = linprog(f,A2,b2,[],[],lb,ub)

%---------------------------------------------------------------
% PROBLEM 1 - PART 3
%---------------------------------------------------------------
% updated constraint matrix A3 with <=25% B-rated bond row
A2 = [-10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -10 -7 -100 1 0 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -10 -8 -6 -110 -107 0 -1.02 1 0 0 0;
     -10 -7 -8 -6 -7 -6 -5 -110 -108 -106 0 0 0 0 -1.03 1 0 0;
     -10 -7 -8 -6 -7 -106 -105 0 0 0 0 0 0 0 0 -1.04 1 0;
     -10 -7 -8 -106 -107 0 0 0 0 0 0 0 0 0 0 0 -1.05 1; 
     -100 -107 -108 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.06;
     324 282 297 278.1 289.8 287.7 -92.9 -110 -104 -101 -107 
     -102 -95.2 0 0 0 0 0];

% updated constraint vector b3
b2 = [-500; -200; -800; -400; -700; -900; 0];

% calculate new objective function value fval3 and dec. var x3
[x3, fval3] = linprog(f,A2,b2,[],[],lb,ub)

%----------------------------------------------------------------
% PROBLEM 2 - PART 1A 
%----------------------------------------------------------------
% Load .csv data from yahoo finance
SPY1 = readmatrix('SPY.csv');    % Jan 14 to Jan 21 monthly data
SPY2 = readmatrix('SPY2.csv');   % Feb 14 to Feb 21 monthly data
GOVT1 = readmatrix('GOVT.csv');  % Jan 14 to Jan 21 monthly data
GOVT2 = readmatrix('GOVT2.csv'); % Feb 14 to Feb 21 monthly data
EEMV1 = readmatrix('EEMV.csv');  % Jan 14 to Jan 21 monthly data
EEMV2 = readmatrix('EEMV2.csv'); % Feb 14 to Feb 21 monthly data

% Return (r_it) = (Month 2 Close - Month 1 Close)/Month 1 Close
SPY_rtn = (SPY2(:, 6) - SPY1(:,6))./SPY1(:, 6);      % mo. return
GOV_rtn = (GOVT2(:, 6) - GOVT1(:,6))./GOVT1(:, 6);   % mo. return
EEM_rtn = (EEMV2(:, 6) - EEMV1(:,6))./EEMV1(:, 6);   % mo. return

% Calculate arithmetic means (rbar_i)
SPY_rbar = mean(SPY_rtn)
GOV_rbar = mean(GOV_rtn)
EEM_rbar = mean(EEM_rtn)

% Calculate geometric Means (mu_i)
SPY_rtn_t = 1;      % initialize
GOV_rtn_t = 1;      % intiialize
EEM_rtn_t = 1;      % initialize

% loop through to calculate total product of returns
for t = 1:85        %T = 85 data points
    SPY_rtn_t = (1 + SPY_rtn(t,:))*SPY_rtn_t;   % cumulative
    GOV_rtn_t = (1 + GOV_rtn(t,:))*GOV_rtn_t;   % cumulative
    EEM_rtn_t = (1 + EEM_rtn(t,:))*EEM_rtn_t;   % cumulative
end

% take the 'T'th root to get final geometric mean (mu_i)
SPY_mu = SPY_rtn_t^(1/85) - 1
GOV_mu = GOV_rtn_t^(1/85) - 1
EEM_mu = EEM_rtn_t^(1/85) - 1

% Calculate standard deviations
SPY_std = std(SPY_rtn)
GOV_std = std(GOV_rtn)
EEM_std = std(EEM_rtn)

% Calculate covariances (SPY: i=1, GOV: i=2, EEM: i=3)
rtn_matrix = [SPY_rtn GOV_rtn EEM_rtn];     % i
format shortE
cov_matrix = cov(rtn_matrix,1)

%-----------------------------------------------------------------
% PROBLEM 2 - PART 1B
%----------------------------------------------------------------
% Generating efficient frontier of three assets: SPY, GOVT, EEMV with
% shorting

R = linspace(0.03, 0.18, 15);    % expected return data points
H = cov_matrix;                  % covariance matrix
f = zeros(3,1);                  % [0 0 0]T

Aeq = [-SPY_mu -GOV_mu -EEM_mu; ones(1,3)]  %Equal. constr. matrix

% initiate vector, matrix to store values from loop:
p_var = zeros(size(R,1),1);     % portfolio variance (obj function)
w = zeros(3, size(R,1));        % asset weights (dec. variable)

% loop through i from 1 to 15
for i = 1:size(R,2)
    beq = [-R(i); 1];                % Equality constraint vector
    % quadratic program
    [w(:,i), p_var(i)] = quadprog(H, f, [], [], Aeq, beq,[],[]);
end

% Results table of: 'R', 'p_var', w_1, w_2, w_3:
T = array2table([R',p_var', w']);
T.Properties.VariableNames(1:5) = {'exp. rtn', 'port. var', 'w_1', 'w_2', 'w_3'}
display(T)

% % Plot efficient frontier
% figure('Name','Efficient Frontier - With Shorting');
% plot(p_var, R, 'b-*');
% title('Figure 1 - Efficient Frontier with SPY, GOVT, EEMV');
% xlabel('portfolio variance');
% ylabel('expected return');

%-----------------------------------------------------------------
% PROBLEM 2 - PART 1C
%----------------------------------------------------------------
% minimum variance portfolio: from Table T, occurs at: 
% exp_rtn = 0.03, var = 0.006, w_1 = 3.9, w_2 = -0.04, w_3 = -2.9

% use feb_rtns for exp rtn calculations
feb_rtns = [SPY_rtn(85,:) GOV_rtn(85,:) EEM_rtn(85,:)];

% equal weight portfolio
w3 = [0.33; 0.33; 0.33];
var_3 = w3' * H * w3
exp_rtn_3 = feb_rtns * w3

% 60% SPY, 30% GOVT, 10% EEMV portfolio
w4 = [0.6; 0.3; 0.1];
var_4 = w4' * H * w4
exp_rtn_4 = feb_rtns * w4

% % Plot efficient frontier
% figure('Name','Efficient Frontier');
% plot(p_var, R, 'b-*')
% hold on
% scatter(0.006, 0.03, 'r*')
% scatter(var_3, exp_rtn_3, 'm*')
% scatter(var_3, exp_rtn_3, 'g*')
% legend('efficient frontier', 'min-var portfolio', 
%        'equal wt portfolio', '60-40-10 portfolio')
% title('Figure 2 - Compare 3 Portfolios');
% xlabel('portfolio variance');
% ylabel('expected return');

%----------------------------------------------------------------
% PROBLEM 2 - PART 2A 
%----------------------------------------------------------------
% Load .csv data from yahoo finance
CME1 = readmatrix('CME.csv');    % Jan 14 to Jan 21 monthly data
CME2 = readmatrix('CME2.csv');   % Feb 14 to Feb 21 monthly data
BR1 = readmatrix('BR.csv');      % Jan 14 to Jan 21 monthly data
BR2 = readmatrix('BR2.csv');     % Feb 14 to Feb 21 monthly data
CBOE1 = readmatrix('CBOE.csv');  % Jan 14 to Jan 21 monthly data
CBOE2 = readmatrix('CBOE2.csv'); % Feb 14 to Feb 21 monthly data
ICE1 = readmatrix('ICE.csv');    % Jan 14 to Jan 21 monthly data
ICE2 = readmatrix('ICE2.csv');   % Feb 14 to Feb 21 monthly data
ACN1 = readmatrix('ACN.csv');    % Jan 14 to Jan 21 monthly data
ACN2 = readmatrix('ACN2.csv');   % Feb 14 to Feb 21 monthly data

% Return (r_it) = (Month 2 Close - Month 1 Close)/Month 1 Close
CME_rtn = (CME2(:, 6) - CME1(:,6))./CME1(:, 6);      % mo. return
BR_rtn  = (BR2(:, 6) - BR1(:,6))./BR1(:, 6);         % mo. return
CBO_rtn = (CBOE2(:, 6) - CBOE1(:,6))./CBOE1(:, 6);   % mo. return
ICE_rtn = (ICE2(:, 6) - ICE1(:,6))./ICE1(:, 6);      % mo. return
ACN_rtn = (ACN2(:, 6) - ACN1(:,6))./ACN1(:, 6);      % mo. return

% Calculate arithmetic means (rbar_i)
CME_rbar = mean(CME_rtn)
BR_rbar  = mean(BR_rtn)
CBO_rbar = mean(CBO_rtn)
ICE_rbar = mean(ICE_rtn)
ACN_rbar = mean(ACN_rtn)

% Calculate geometric Means (mu_i)
CME_rtn_t = 1;      % initialize
BR_rtn_t = 1;      % intiialize
CBO_rtn_t = 1;      % initialize
ICE_rtn_t = 1;      % intiialize
ACN_rtn_t = 1;      % initialize

% loop through to calculate total product of returns
for t = 1:85        %T = 85 data points
    CME_rtn_t = (1 + CME_rtn(t,:))*CME_rtn_t;   % cumulative
    BR_rtn_t = (1 + BR_rtn(t,:))*BR_rtn_t;      % cumulative
    CBO_rtn_t = (1 + CBO_rtn(t,:))*CBO_rtn_t;   % cumulative
    ICE_rtn_t = (1 + ICE_rtn(t,:))*ICE_rtn_t;   % cumulative
    ACN_rtn_t = (1 + ACN_rtn(t,:))*ACN_rtn_t;   % cumulative
end

% take the 'T'th root to get final geometric mean (mu_i)
CME_mu = CME_rtn_t^(1/85) - 1
BR_mu = BR_rtn_t^(1/85) - 1
CBO_mu = CBO_rtn_t^(1/85) - 1
ICE_mu = ICE_rtn_t^(1/85) - 1
ACN_mu = ACN_rtn_t^(1/85) - 1

% Calculate standard deviations
CME_std = std(CME_rtn)
BR_std = std(BR_rtn)
CBO_std = std(CBO_rtn)
ICE_std = std(ICE_rtn)
ACN_std = std(ACN_rtn)

% Calculate covariances 
% (SPY=1, GOV=2, EEM=3,CME=4, BR=5, CBO=6, ICE=7, ACN=8)
rtn_matrix2 = [SPY_rtn GOV_rtn EEM_rtn 
               CME_rtn BR_rtn CBO_rtn ICE_rtn ACN_rtn];
format shortE
cov_matrix2 = cov(rtn_matrix2,1)

%-----------------------------------------------------------------
% PROBLEM 2 - PART 2B
%----------------------------------------------------------------
% Generating efficient frontier using all 8 assets:

R2 = linspace(0.03, 0.33, 10);     % expected return data points
H2 = cov_matrix2;                  % covariance matrix
f2 = zeros(8,1);                    % [0 0 0 0 0 0 0 0]T
 
Aeq2 = [SPY_mu GOV_mu EEM_mu CME_mu BR_mu CBO_mu ICE_mu ACN_mu; 
       ones(1,8)];  %Equality constraint matrix

% initiate vector, matrix to store values from loop:
p_var2 = zeros(size(R2,1),1);  % portfolio variance (obj function)
w2 = zeros(8, size(R2,1));     % asset weights (dec. variable)

% loop through i from 1 to 15
for i = 1:size(R2,2)
    beq2 = [R2(i); 1];          % Equality constraint vector
    % quadratic program
    [w2(:,i), p_var2(i)] = quadprog(H2, f2, [], [], Aeq2, beq2,[],[]);
end

% Results table of: 'R', 'p_var', w_i
T2 = array2table([R2',p_var2', w2']);
T2.Properties.VariableNames(1:10) = {'exp return', 'port. var','w_1', 'w_2', 'w_3', 'w_4', 'w_5', 'w_6', 'w_7', 'w_8'}
display(T2)
% 
% % Plot efficient frontier
% figure('Name','Efficient Frontier - With Shorting');
% plot(p_var2, R2, 'b-*');
% hold on
% plot(p_var, R, 'r--o');
% legend('8-asset portfolio', '3-asset portfolio')
% title('Figure 3 - Efficient Frontier with 8 assets');
% xlabel('portfolio variance');
% ylabel('expected return');
