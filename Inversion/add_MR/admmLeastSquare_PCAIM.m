% Inversion taking into account non-negativity constrain as given by ADMM
% (Alternating Direction Method of Multipliers) .
% PROBLEM TO FIX: does not work when Cm cannot be inverted (is singular)
% The function to minimize is:
% f(m,F) = (m-m0)' * Cm_inv * (m-m0)
%        + (G*m-d)'* Cd_inv * (G*m-d)
%        + rho/2*(F-m+lambda)' * (F-m+lambda)
% with constrain that m(i) >=0 for all i
% This is solved iteratively
% (1)  find m_(1+k) E argmin f(M,F_k)
%      m_k = inv (Gt * Cd_inv * G + Cm_inv + rho/2 * Id) 
%          * (Gt * Cd_inv * dt + Cm_inv * m0 + rho/2 * Gt * (F+lambda))
% (2)  find F_(1+k) E argmin f(M^(1+k),F)
%      F_k(i) = max(0,m_(k+1)(i)-lambda_k(i))
% (3)  lambda_(k+1) = lamda_k + (f_(k+1)- m_(k+1))

% When rho increase, convergence is faster rho ~ [1-100]
% n_iter has to be large. Typically, rho = 100 and n_inter = 2000 works

% Mathilde Radiguet
% Created: Feb. 2015
% Modified: Oct. 2020: when ramp are inverted, no positivity imposed on
% ramp parameters
%%
function m1 = admmLeastSquare_PCAIM(G,Cd,Cm,m0,d,rho,n_iter,n_ramp)
% ititial values (lambda_0 and F_0 are fixed to zero, this could be
% modified)
lambda_0 = zeros(length(m0),1);
F_0 = zeros(length(m0),1);
Sf_c = zeros(1,n_iter);
Sdata_c = zeros(1,n_iter);
Sp_c = zeros(1,n_iter);

% initialisation 
Gt = G';
Cd_inv = inv(Cd);
Cm_inv = inv(Cm);
lambda_k = lambda_0;
F_k = F_0;
X_inv = Gt*Cd_inv*G+Cm_inv+rho/2*eye(length(m0));
Y = Gt*Cd_inv*d+Cm_inv*m0';
X = inv(X_inv);
XY = X*Y;
XGt = X ;
GCmGt = G * Cm * Gt;

%%%%%%%%%%%%%%%%%  ADMM iterations
for c = 1:n_iter
    % step (1) 
    m_k = XY + rho/2 * XGt * (F_k + lambda_k);
    % test
%    m_k = Cm * Gt *inv(Cd + G * Cm*Gt + inv(Gt*Cd_inv)*rho/2*Cm*Gt);
    
    % step (2)
    % Positivity for the slip parameters
    F_k(1:end-n_ramp) = max(0,m_k(1:end-n_ramp)-lambda_k(1:end-n_ramp));
    % No constrain on the ramp parameters
    F_k(end-n_ramp:end) = m_k(end-n_ramp:end)-lambda_k(end-n_ramp:end);    
    % step (3) compute new lambda
    lambda_k = lambda_k + (F_k - m_k);
    
end
m1=m_k;
  