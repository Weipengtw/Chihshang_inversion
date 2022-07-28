function [L,L_var,R_diag,R_rest]= invert_components_tarantola(U,G_projected_dense,fault_position,n_comp,options)

%INVERT_COMPONENTS   Set up and execute the inversion for slip at depth
%   L=INVERT_COMPONENTS(U,G_PROJECTED_DENSE,FAULT_POSITION,N_COMP,OPTIONS) finds
%   a possibly constrained, regularized least-squares solution to inverting
%   U for dislocation at depth with Greens functions G_projected_dense.
%   Smoothing is introduced via regularisation in the Cm covariance matrix
%   and a regularisation length lambda, as in Radiguet et al. 2011 GJI
%
%   Example:
%   PCAIM_driver
%
%   See also PCAIM_driver.

%   Written by Mathilde Radiguet
%   Addapted from invert_components.m by
%   Andrew Kositsky and Hugo Perfettini

% cost function = (G * m - d)^t * C_d ^{-1} * (G * m - d)
%               + (m - m_0)^ * C_m^{-1} * (m - m0)

% X=[nstat,nmeas]
% G=[nstat,n_patches]
% G-1=[n_patches,nstat]
% slip=[n_patches,nmeas]
% V=[nmeas,nmeas]
% U=[nstat,nstat]
% S=[nstat,nmeas]
% X=G*slip=U*S*V'
% slip=(G^-1*U)*S*V';
% L=G-1*U=[n_patches,nstat][nstat,nstat]=[n_patches,nstat]
% $L = [G^{-1} U_1, G^{-1} U_2, \cdots, G^{-1} U_r]$
% strike_vector_index_trg=16:18;
% updip_vector_index_trg=19:21;
% strike_vector_index_rect=20:22;
% updip_vector_index_rect=23:25;

% cost function_i = (G * L_i - U_i)^t * C_d ^{-1} * (G * L_i - U_i)
%               + (L_i - m_i0)^ * C_m^{-1} * (L_i - m_i0)

% Invert solution:
% in model space
% L = m0 + (Gt * Cd^{-1} * G + C_m^{-1})^{-1} * Gt * Cd^{-1} * (U - G * m0)
% in data space (no need to compute Cm_inv)
% L = m0 + (Cm * Gt / (G * Cm * Gt + Cd) * (U'-G * m0))

% Define G by concatonating together the Greens functions for the
% individual datasets
%G = extract_from_cell(G_projected_dense);
G = G_projected_dense{1};
% calculate the number of time series and dislocations * dislocaiton_directions
[n_tseries,n_dislocations_dirs]=size(G);

%%% THE FOLLOWING ASSUMPTIONS ARE ONLY TRUE FOR 2-COMPONENT DISLOCATIONS AT
%%% DEPTH, NOT THE GENERAL CASE
% Ensure that the number of dislocations * dislocation_dirs is even, then
% set the number of patches
if mod(n_dislocations_dirs,2)~=0,
    error('The number of Columns of G should be even');
end
n_patches=n_dislocations_dirs/2;

%Determine the Number of Options
n_options=size(options,2);

%Default Values for all non-mandatory options
V_sign = 1;
U_var = speye(numel(U));
Em0=ones(1,n_patches);
fixed_rake_flag=0;
L_var_flag=0;

%Load Options
kk = 1;
while kk <= n_options,
    switch options{kk}
        case 'Vsign' % the sign of V in the first component (for positivity)
            V_sign = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the V_sign value
        case 'FixedRake'  % Force fixed rake on the fault
            fixed_rake_flag=1;
            rake=options{1,kk+1};
            if ~isnumeric(rake),
                error('Rake angle was not given!');
            end
            norm_flag=1;%When using fixed rake, the normalization of slip is OBLIGATORY
            kk = kk+1; %Extra increment so we skip the rake value
        case 'Uvar'
            U_var = options{1,kk+1};
            kk = kk+1;
        case 'ComputeLvar'
            L_var_flag = 1;
        case 'WhichSmoothing'
            which_smoothing = options{1,kk+1};
            kk = kk+1;
        case 'LambdaStrike'
            lambda_strike =  options{1,kk+1};
            kk = kk+1;
        case 'LambdaDip'
            lambda_dip =  options{1,kk+1};
            kk = kk+1;
        case 'Lambda0'
            lambda0 =  options{1,kk+1};
            kk = kk+1;
        case 'WhichSpace'
            which_space = options{1,kk+1};
            kk = kk+1;
        case 'Em0'
            Em0 = options{1,kk+1};
            kk = kk+1;
        otherwise,
            if ischar(options{1,kk}),
                sub_str = options{1,kk};
            elseif isnumeric(options{1,kk})
                sub_str = num2str(options{1,kk});
            end
            error([sub_str,' is not a valid option.']);
    end
    kk = kk+1; %Increment to next value
end

if ~isempty(U_var),
    disp('Uvar ponderation');
    Cd = U_var;
else
    disp('No Uvar ponderation');
end


% Project the Greens function onto rake vectors if fixed_rake_flag ==1.
% Otherwise, assign G_proj = G.

if fixed_rake_flag ==1
    G_proj = project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag,V_sign);
    G_proj_blk = blkdiag_rep(G_proj,n_comp);
else
    G_proj_blk = blkdiag_rep(G,n_comp);
end

% compute Cm
% Set matrix of block diagonal repetition of Cm
if fixed_rake_flag ==1 % only half of the matrix is non-zero
    Cm_non_null = compute_covariance_matrix_model(fault_position,n_patches,Em0,lambda_strike,lambda0,which_smoothing);
    Cm_null = zeros(n_patches,n_patches);
    Cm = blkdiag(Cm_non_null,Cm_null);
    Cm_blk = blkdiag_rep(Cm,n_comp);
else % all the matrix is non-zero
    Cm_non_null_strike = compute_covariance_matrix_model(fault_position,n_patches,Em0,lambda_strike,lambda0,which_smoothing);
    Cm_non_null_dip = compute_covariance_matrix_model(fault_position,n_patches,Em0,lambda_dip,lambda0,which_smoothing);
    Cm = blkdiag(Cm_non_null_strike,Cm_non_null_dip);
    Cm_blk = blkdiag_rep(Cm,n_comp);
end

% Set initial model concatenated
m0 = zeros(1,n_patches*2);
m0_blk = repmat(m0,1,n_comp);

% data vector
d = U(:);

% Set the initial design matrix equal to a block-diagonal repetition of the
% projected Greens function
G_proj_blk_t = G_proj_blk';
switch which_space
    case 'model_space' % invert in the model space
        Cm_blk_inv = inv(Cm_blk);
        Cd_inv = inv(Cd);
        A = (G_proj_blk_t * Cd_inv * G_proj_blk_t + Cm_blk_inv);
        % compute final model
        m1 = m0 + A\G_proj_blk_t * Cd_inv * (d.' - G_proj_blk * m0);
        % compute final covariance matrix
        if L_var_flag,
            Cm_1 = inv(G_proj_blk_t * Cd_inv * G_proj_blk + Cm_blk_inv);
        else
            disp('L_var is not computed. Output gives L_var=0');
        end
    case 'data_space' % invert in the data space (no need to compute Cm_inv)
        % compute final model
        B = inv(G_proj_blk * Cm_blk * G_proj_blk_t + Cd);
        m1 = m0_blk' + (Cm_blk * G_proj_blk_t / (G_proj_blk * Cm_blk * G_proj_blk_t + Cd) * (d - G_proj_blk * m0_blk'));
        % compute final covariance matrix
        if L_var_flag,
            Cm_1 = Cm_blk - Cm_blk * G_proj_blk_t * B * G_proj_blk * Cm_blk;
        else
            disp('L_var is not computed. Output gives L_var=0');
        end
    otherwise
        error([which_space 'is not a valid option for regularLeastSquare. Choices are model_space and data_space'])
end

% Compute resolution matrix
R = Cm_blk * G_proj_blk_t / ( G_proj_blk * Cm_blk * G_proj_blk_t + Cd) * G_proj_blk;
R_diag = reshape(diag(R),n_patches*2,n_comp);
R_rest = reshape(sum(R'),n_patches*2,n_comp);

% recompose slip vector in case of fixed rake
if fixed_rake_flag==1,
    m1 = recompose_slip_vectors(rake,n_patches,n_comp,m1,V_sign,norm_flag);
end

% Reshape the slip at depth vector so that G*L = U is approximately true.
L = reshape(m1,n_patches*2,n_comp);
if L_var_flag,
    L_var = Cm_1; % CORRECT NOW
else
    L_var=0;
end
disp('Inversion Finished');

function G_proj = project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag,V_sign)
%PROJECT_GREEN_FUNCTION   Project Greens functions onto a fixed rake vector
%   G_PROJ = PROJECT_GREEN_FUNCTION(G,N_PATCHES,FIXED_RAKE_FLAG,RAKE,NORM_FLAG)
%   projects the Greens functions G onto rake directions RAKE if
%   FIXED_RAKE_FLAG == 1. Otherwise, G_PROJ =G.
%
%   Example:
%
%   PCAIM_driver
%
%   See also INVERT_COMPONENTS, PCAIM_DRIVER.

%   By Hugo Perfettini
%   Commented and revised by Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/12/03  $

% Divide up the original Greens functions into strike-slip (G1) and
% dip-slip (G2) components.
G1=G(:,1:n_patches);
G2=G(:,n_patches+1:end);

% If the fixed_rake_flag == 1, then continue.
if fixed_rake_flag==1,
    %%%% HUGO'S NEW VERSION (edited by APK post-hoc)
    %Preallocate the projected Greens function matrix G_proj
    G_proj = zeros(size(G));
    % Calculate the number of observation points on the surface
    n_points=size(G1,1);
    % For each dislocation source
    pure_strike_slip_index = find(mod(rake,180)==0);
    pure_dip_slip_index = find(mod(rake,180)== 90);
    non_pure_index = setdiff(setdiff((1:n_patches),pure_dip_slip_index),pure_strike_slip_index);
    
    ds_pure_index = [pure_dip_slip_index,pure_dip_slip_index+n_patches];
    ss_pure_index = [pure_strike_slip_index,pure_strike_slip_index+n_patches];
    ss_and_ds_non_pure_index = [non_pure_index,non_pure_index+n_patches];
    
    %Correct direction factor so that positive slip corresponds to the
    %   proper direction of slip
    if numel(ds_pure_index) > 0
        G_proj(:,ds_pure_index) = ...
            V_sign*[zeros(size(G1(:,pure_dip_slip_index)))   ,...
            repmat(sign(sind(rake(pure_dip_slip_index)))',size(G2,1),1).*...
            G2(:,pure_dip_slip_index)];
    end
    if numel(ss_pure_index) > 0
        G_proj(:,ss_pure_index) = ...
            V_sign*[repmat(sign(cosd(rake(pure_strike_slip_index)))',size(G1,1),1).*G1(:,pure_strike_slip_index),...
            zeros(size(G2(:,pure_strike_slip_index)))];
    end
    G_proj(:,ss_and_ds_non_pure_index) =...
        V_sign*[repmat(sign(cosd(rake(non_pure_index)))',size(G1,1),1) .*...
        G1(:,non_pure_index)+...
        repmat(tand(rake(non_pure_index))'.*sign(cosd(rake(non_pure_index)))'...
        , size(G2,1),1)...
        .* G2(:,non_pure_index),...
        zeros(size(G1(:,non_pure_index)))];
    
    if norm_flag==1,
        G_proj(:,ss_and_ds_non_pure_index)=G_proj(:,ss_and_ds_non_pure_index)./ ...
            ...            repmat(sqrt(1+tand([rake(non_pure_index);rake(non_pure_index)]).^2)',n_points,1)*sqrt(2);%Normalization of Green Functions
            repmat(sqrt(1+(tand([rake(non_pure_index);rake(non_pure_index)]).^2))',n_points,1);%Normalization of Green Functions
    end
else
    % Otherwise, just assign G_proj = [G1,G2] = G
    G_proj=[G1,G2];
end

function B = blkdiag_rep(A,n)
eval(['B = blkdiag(' repmat('A, ',1,n-1),'A);']);