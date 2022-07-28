function Cm = compute_covariance_matrix_model(fault_model,n_patches,Em0,lambda,lambda0,which_smoothing)
% Written by Mathilde Radiguet
% Updated 2020/10/27 to allow variable lambda_0 and allow variability in
% sigma_m (Em_sq(i,j) = Em(i) * Em(j))
%COMPUTE_COVARIANCE_MATRIX_MODEL computes the covariance matrix for model parameters
% The smoothing is introduced at this step
% Model covariance matric contains
% Cm(i,j) = Em(i)^2 * exp( - 0.5 * d(i,j)^2 / lambda^2) (gaussian)
% Cm(i,j) = Em(i)^2 * exp(- d(i,j) / lambda) (exponential)
% where
% d(i,j) = euclidian distance between barycenter of patches i and j
% Em(i) = weighting factor, related to the smoothing
% (Valette, 05/11/2009), and see Rediguet et al., 2011, GJI
% Em(i) = Em0(i).*lambda0 ./lambda
% lambda is the caracteristic length of regularisation
% lambda0 is a scaling factor 
% Em0(i) is the initial standard devation on model parameters
% Should be the sam everywhere, unless one wants to prevent slip on a
% specific fault patch.

% Em0(1,n_patches)
% Em(1,npatches)
% Cm(npatches,npatches)

% Normalize standard deviation on model depending on the smoothing (see Valette, 05/11/2009)
%Em=Em0.*lambda0'./lambda;
Em=Em0.*lambda0'./lambda;
% initialize
Cm=zeros(n_patches,n_patches);
% compute distance from fault model
d=zeros(n_patches,n_patches);
for i = 1:length(fault_model(:,1))
    for j = 1:length(fault_model(:,1))
        d(i,j) = sqrt((fault_model(i,1)-fault_model(j,1))^2+(fault_model(i,2)-fault_model(j,2))^2+(fault_model(i,3)-fault_model(j,3))^2);
        Em_sq(i,j) = Em(i) * Em(j);
    end
end
%fact=sum(sum(ss));
%ss=ss/fact;
switch which_smoothing
    case 'gaussian' 
        % use bsxfun to multiply each line of Cm by Em
        Cm = bsxfun(@times,exp(-0.5*d.^2/lambda^2),Em_sq);
        %Cm = bsxfun(@times,exp(-0.5*d.^2/lambda^2),Em.^2);
    case 'exponential'
        Cm = bsxfun(@times,exp(-d/lambda),Em_sq);
        %Cm = bsxfun(@times,exp(-d/lambda),Em.^2);
    otherwise
        error([which_smoothing 'is not a valid option for smoothing choice. Choices are expotential_smoothing and gaussian_smoothing'])
end

   