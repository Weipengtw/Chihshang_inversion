function [L,L_var]=invert_components_combined_NEW(U,G_projected_dense,gamma,Lap,n_comp,options)
%INVERT_COMPONENTS   Set up and execute the inversion for slip at depth
%   L=INVERT_COMPONENTS(U,G_PROJECTED_DENSE,GAMMA,LAP,N_COMP,OPTIONS) finds
%   a possibly constrained, regularized least-squares solution to inverting
%   U for dislocation at depth with Greens functions G_projected_dense.
%   Sparse constraints can be added in options. GAMMA is the smoothing of
%   the Laplacian LAP, and N_COMP is the number of components in the model.
%
%   Example:
%   PCAIM_driver
%
%   See also PCAIM_driver.

%   Rewritten to only invert once.

%   By Andrew Kositsky and Hugo Perfettini
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.1.0.0 $  $Date: 2010/06/30  $

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

%   Written by Hugo Perfettini
%   Commented and revised by Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/11/30  $

% Define G by concatonating together the Greens functions for the
% individual datasets
G=extract_from_cell(G_projected_dense);
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
V_sign = 1;positivity_flag=0;fixed_rake_flag=0;pinv_flag=0;norm_flag=0;lap_flag=1;
fnnls_flag=1;lsqnonneg_flag=0;direction_factor=1;rake=NaN;lsqlin_options_flag=0;linprog_options_flag=0;
epsilon = 10^(-7); U_var = speye(numel(U));

if gamma==0,lap_flag=0;disp('No smoothing applied');end

%Load Options
kk = 1;
while kk <= n_options,
    switch options{kk}
        case 'Vsign' % the sign of V in the first component (for positivity)
            V_sign = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the V_sign value
        case 'Positivity' % Force positive slip at depth with one component
            positivity_flag=1;
        case 'lsqnonneg' % Use lsqnonneg instead of lsqlin (for positivity)
            lsqnonneg_flag=1;
        case 'FixedRake'  % Force fixed rake on the fault
            fixed_rake_flag=1;
            rake=options{1,kk+1};
            if ~isnumeric(rake),
                error('Rake angle was not given!');
            end
            norm_flag=1;%When using fixed rake, the normalization of slip is OBLIGATORY
            kk = kk+1; %Extra increment so we skip the rake value
        case 'PseudoInverse'
            pinv_flag=1;
        case 'NormalizeGreenFunctions'
            norm_flag=1;
        case 'SparseConstraint'
            SparseConstraint=options{1,kk+1};
            if(~isnumeric(SparseConstraint) &&...
                    size(SparseConstraint,2) == n_comp*n_patches*2)
                error(['SparseConstraint must be a N x '...
                    num2str(n_comp*n_patches*2) 'matrix']);
            end
            kk = kk+1; %Extra increment so we skip the sparse constrant matrix
        case 'SparseWeight'
            SparseWeight = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the sparse weight vector
        case 'Sparse_d'
            Sparse_d = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the sparse data vector
        case 'NoSmoothing'
            lap_flag=0;%Flip the flag for no smoothing on the fault surface
        case 'lsqlin_options'
            lsqlin_options_flag=1;
            lsqlin_options=options{1,kk+1};
            kk = kk+1;
        case 'linprog_options'
            linprog_options_flag=1;
            linprog_options=options{1,kk+1};
            kk = kk+1;
        case 'epsilon'
            epsilon=options{1,kk+1};
            kk = kk+1;
        case 'Uvar'
            U_var = options{1,kk+1};
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
    U_var_inv = U_var^(-1);
    U_var_inv_chol= chol(U_var_inv);
else
    disp('No Uvar ponderation');
    U_var_inv_chol=eye(size(U));
end

G_proj=project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag,V_sign);
G_proj_combined=blkdiag_rep(G_proj,n_comp);

d_tot=U_var_inv_chol*U(:);
A_tot=U_var_inv_chol*G_proj_combined;

% Add the Laplacian to the design matrix and zeros to the data vector if
% using the Laplacian for regularization.

if lap_flag==1,
    d_combined_tmp=[];
    A_combined_tmp=[];
    if n_comp==2,
        for ii=1:n_comp,
            d_combined_tmp=[d_combined_tmp;zeros(size(Lap{ii},1),1)];
            A_combined_tmp=blkdiag(A_combined_tmp,gamma*Lap{ii});
        end
    else
        d = [d_combined_tmp;zeros(size(Lap{1},1)*n_comp,1)];
        A = [A_combined_tmp;blkdiag_rep(Lap{1},n_comp)*gamma];
    end
end

save A_COMBINED_NEW.mat  A_tot U_var_inv_chol  G_proj_combined G_proj;

d_combined=[d_tot;d_combined_tmp];
A_combined=[A_tot;A_combined_tmp];

d=d_combined;
A=A_combined;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there exists a sparse data-based constraint (defined in options), then
% add this to the augmented design matrix and the data vector

if exist('SparseConstraint','var')
    n_sparse_datasets=numel(SparseConstraint);
    for kk=1:n_sparse_datasets,
        if SparseWeight(kk)>0,
            SparseConstraint_tmp = project_Green_function(SparseConstraint{kk},...
                n_patches,fixed_rake_flag,rake,norm_flag,V_sign);
            A = [A;SparseConstraint_tmp*SparseWeight(kk)];
            d = [d;Sparse_d{kk}*SparseWeight(kk)];
        end
    end
end
%%% WARNING
%%% WARNING
%%% WARNING
%%% THIS SECTION ONLY WORKS FOR TWO-COMPONENT DISLOCATIONS (STRIKE + DIP)
%%% AND IS NOT GENERAL TO ALL TYPES OF DISLOCATION.
%%%

% Create indexes for the strike-slip and dip-slip patches,
%%% CURRENTLY THESE ARE UNUSED
%%%repeated_slip_patch_offsets =repmat((0:n_comp-1)*n_patches*2,n_patches,1);
%%%ss_patch_indexes = repmat(1:n_patches,n_comp,1)' + repeated_slip_patch_offsets;
%%%ds_patch_indexes = n_patches + ss_patch_indexes;0

% Define the inversion options set above
inversion_opt = ...
    {positivity_flag,direction_factor,fnnls_flag,...
    pinv_flag,fixed_rake_flag,norm_flag,lsqnonneg_flag,rake,V_sign};
%add lsqlin options

if lsqlin_options_flag,
    n_inversion_opt=numel(inversion_opt);
    [~,ind_lsqlin]=find_string_in_cell(options,'lsqlin');
    inversion_opt{n_inversion_opt+1}=options{ind_lsqlin+1};
end
%add linprog options
if linprog_options_flag,
    n_inversion_opt=numel(inversion_opt);
    [~,ind_linprog]=find_string_in_cell(options,'linprog');
    inversion_opt{n_inversion_opt+1}=options{ind_linprog+1};
end
% Perform the inversion step within inversion_type
s=inversion_type(d,A,n_patches,n_comp,inversion_opt);

L     = reshape(s,n_patches*2,n_comp);

% s_var = (A'*A)^(-1);
%L_var = reshape(s_var,n_patches*2,n_comp);
% L_var=s_var;
L_var=[];

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

% Divide up the original Greens functions into strike-slip (G1_proj) and
% dip-slip (G2_proj) components.
G1_proj=G(:,1:n_patches);
G2_proj=G(:,n_patches+1:end);

% If the fixed_rake_flag == 1, then continue.
if fixed_rake_flag==1,
    
    %%%% HUGO'S NEW VERSION (edited by APK post-hoc)
    %Preallocate the projected Greens function matrix G_proj
    G_proj = zeros(size(G));
    % Calculate the number of observation points on the surface
    n_points=size(G1_proj,1);
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
            V_sign*[zeros(size(G1_proj(:,pure_dip_slip_index)))   ,...
            repmat(sign(sind(rake(pure_dip_slip_index)))',size(G2_proj,1),1).*...
            G2_proj(:,pure_dip_slip_index)];
    end
    if numel(ss_pure_index) > 0
        G_proj(:,ss_pure_index) = ...
            V_sign*[repmat(sign(cosd(rake(pure_strike_slip_index)))',size(G1_proj,1),1).*G1_proj(:,pure_strike_slip_index),...
            zeros(size(G2_proj(:,pure_strike_slip_index)))];
    end
    G_proj(:,ss_and_ds_non_pure_index) =...
        V_sign*[repmat(sign(cosd(rake(non_pure_index)))',size(G1_proj,1),1) .*...
        G1_proj(:,non_pure_index)+...
        repmat(tand(rake(non_pure_index))'.*sign(cosd(rake(non_pure_index)))'...
        , size(G2_proj,1),1)...
        .* G2_proj(:,non_pure_index),...
        zeros(size(G1_proj(:,non_pure_index)))];
    
    if norm_flag==1,
        G_proj(:,ss_and_ds_non_pure_index)=G_proj(:,ss_and_ds_non_pure_index)./ ...
            ...            repmat(sqrt(1+tand([rake(non_pure_index);rake(non_pure_index)]).^2)',n_points,1)*sqrt(2);%Normalization of Green Functions
            repmat(sqrt(1+(tand([rake(non_pure_index);rake(non_pure_index)]).^2))',n_points,1);%Normalization of Green Functions
    end
else
    % Otherwise, just assign G_proj = [G1_proj,G2_proj] = G
    G_proj=[G1_proj,G2_proj];
end

function B = blkdiag_rep(A,n)
eval(['B = blkdiag(' repmat('A, ',1,n-1),'A);']);