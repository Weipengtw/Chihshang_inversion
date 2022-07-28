function [L,L_var,SARramp]=invert_components_SAR_NEW(U,G_projected_dense,gamma,Lap,n_comp,options,...
    all_position,n_dense)
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
if ~isempty(G_projected_dense)	%coseismic using only SAR data
    G = extract_from_cell(G_projected_dense);
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
else
    n_patches=size(Lap,1)/2;
    G=[];
end

%Determine the Number of Options
n_options=size(options,2);

%Default Values for all non-mandatory options
V_sign = 1;positivity_flag=0;fixed_rake_flag=0;pinv_flag=0;norm_flag=0;lap_flag=1;
fnnls_flag=1;lsqnonneg_flag=0;direction_factor=1;rake=NaN;lsqlin_options_flag=0;
epsilon = 10^(-7); U_var=[];SparseFringe_flag=0;MinSumSlip_flag=0;
if gamma==0,lap_flag=0;disp('No smoothing applied');end
[iflag,index]=find_string_in_cell(options,'lsqlin_options');
if iflag==1
    lsqlin_options_flag=1;
    lsqlin_options=options{1,index+1};
end
[iflag,index]=find_string_in_cell(options,'Positivity');
if iflag==1
    positivity_flag=1;
    if lsqlin_options_flag==1
        [iflag,index]=find_string_in_cell(lsqlin_options,'lb');
        if iflag==1,
            lsqlin_options{index+1}=ones(size(lsqlin_options{index+1}))*0;
        else
            lsqlin_options={'lb',ones(1,n_patches*2)*0};
        end
    end
end
%Load Options
kk = 1;
while kk <= n_options,
    switch options{kk}
        case 'Vsign' % the sign of V in the first component (for positivity)
            V_sign = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the V_sign value
        case 'Positivity' % Force positive slip at depth with one component
        case 'lsqnonneg' % Use lsqnonneg instead of lsqlin (for positivity)
            lsqnonneg_flag=1;
        case 'FixedRake'  % Force fixed rake on the fault
            fixed_rake_flag=1;
            rake=options{1,kk+1};
            if ~isnumeric(rake),
                error('Rake angle was not given!');
            end
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
            kk = kk+1;
        case 'SparseWeight'
            SparseWeight = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the sparse weight vector
        case 'Sparse_d'
            Sparse_d = options{1,kk+1};
            kk = kk+1; %Extra increment so we skip the sparse data vector
        case 'NoSmoothing'
            lap_flag=0;%Flip the flag for no smoothing on the fault surface
        case 'lsqlin_options'
            kk = kk+1;
        case 'epsilon'
            epsilon=options{1,kk+1};
            kk = kk+1;
        case 'Uvar'
            U_var = options{1,kk+1};
            kk = kk+1;
        case 'SparseFringe'
            SparseFringe_flag=1;
            SparseFringe=options{1,kk+1};
            if(~isnumeric(SparseFringe) || numel(SparseFringe) ~= 6)
                error('SparseFringe must be a vector of 6 components');
            end
            kk=kk+1;
            if lsqlin_options_flag==1,
                [iflag,index]=find_string_in_cell(lsqlin_options,'lb');
                lb=lsqlin_options{index+1};
                if size(lb,2)==1, lb=lb';end
                if iflag ,lsqlin_options{index+1}=[lb,SparseFringe(1:2:end)];end
                [iflag,index]=find_string_in_cell(lsqlin_options,'ub');
                ub=lsqlin_options{index+1};
                if size(ub,2)==1, ub=ub';end
                if iflag ,lsqlin_options{index+1}=[ub,SparseFringe(2:2:end)];end
            else
                lsqlin_options_flag=1;
                lsqlin_options={'lb',[ones(1,n_patches*2)*-inf,SparseFringe(1:2:end)],'ub',[ones(1,n_patches*2)*inf,SparseFringe(2:2:end)]};
            end
        case 'MinSumSlip'
            
            MinSumSlip_flag=1;
            MinSumSlip=options{1,kk+1};
            kk=kk+1;
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
    %     U_var_inv_chol=eye(size(U));
    U_var_inv_chol=speye(numel(U));
end

% Project the Greens function onto rake vectors if fixed_rake_flag ==1.
% Otherwise, assign G_proj = G.
% G_proj = project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag);
if ~isempty(G)
    G_proj = project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag,V_sign);
    
    % Set the data vector equal to cholsky decomposition of the inverse of the
    %   covariance multiplied by [U(:,1); U(:,2); ... U(:,n_comp)]
    d=U_var_inv_chol*U(:);%
    
    % Set the initial design matrix equal to a block-diagonal repetition of the
    % projected Greens function
    A = U_var_inv_chol*blkdiag_rep(G_proj,n_comp);
end
%Add the Laplacian to the design matrix and zeros to the data vector if
%using the Laplcian for regularization.
if lap_flag==1,
    d = [d;zeros(size(Lap,1)*n_comp,1)];
    A = [A;blkdiag_rep(Lap,n_comp)*gamma];
end

%Add the constraint on the sum of the slip
if MinSumSlip_flag==1
    % A=[A;MinSumSlip(1)*ones(1,size(A,2)/2),zeros(1,size(A,2)/2)];
    A=[A;MinSumSlip(1)*ones(1,size(A,2)/2),ones(1,size(A,2)/2)];
    d=[d;MinSumSlip(1)*MinSumSlip(2)];
end
%% Compactness a tester ne marche pas très bien
% load fault_model
% Ac=zeros(n_patches,n_patches*2);
% for ii=1:n_patches
%     for jj=1:n_patches
%         vec=[fault_model(ii,1)-fault_model(jj,1),fault_model(ii,2)-fault_model(jj,2),...
%             fault_model(ii,3)-fault_model(jj,3)];
%         Ac(ii,jj)=norm(vec)+1;
%         Ac(ii,jj+n_patches)=norm(vec)+1;
%     end
% end
% Ac=Ac/max(max(Ac));
% beta=0.0001;
% A=[A;Ac*beta];
% d=[d;beta*ones(n_patches,1)];
kkramp=0;
if SparseFringe_flag==1
    kkramp=numel(SparseConstraint)*3;
    A = [A ,zeros(size(A,1),kkramp)];
end
% If there exists a sparse data-based constraint (defined in options), then
% add this to the augmented design matrix and the data vector

if exist('SparseConstraint','var'),
    n_sparse_datasets=numel(SparseConstraint);
    kir=1;
    for kk=1:n_sparse_datasets,
        if SparseWeight(kk)>0,
            SparseConstraint_tmp = project_Green_function(SparseConstraint{kk},...
                n_patches,fixed_rake_flag,rake,norm_flag,V_sign);
            if SparseFringe_flag==1
                long_sparse=all_position{kk+n_dense}(:,1);lat_sparse=all_position{kk+n_dense}(:,2);
                xy_sparse=llh2localxy([lat_sparse,long_sparse]',[mean(lat_sparse),mean(long_sparse)]);
                C=zeros(size(SparseConstraint_tmp,1),kkramp);
                C(:,kir:kir+2)=[xy_sparse(:,1)*SparseWeight(kk),xy_sparse(:,2)*SparseWeight(kk),ones(size(xy_sparse(:,2)))*SparseWeight(kk)];
                A = [A;SparseConstraint_tmp*SparseWeight(kk),C];
                d = [d;Sparse_d{kk}*SparseWeight(kk)];
                kir=kir+3;
            else
                A = [A;SparseConstraint_tmp*SparseWeight(kk)];
                d = [d;Sparse_d{kk}*SparseWeight(kk)];
            end
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
%%%ds_patch_indexes = n_patches + ss_patch_indexes;
% Define the inversion options set above
inversion_opt = ...
    {positivity_flag,direction_factor,fnnls_flag,...
    pinv_flag,fixed_rake_flag,norm_flag,lsqnonneg_flag,rake,V_sign};
%add lsqlin options
n_inversion_opt=numel(inversion_opt);
if lsqlin_options_flag,inversion_opt{n_inversion_opt+1}=lsqlin_options;end
% Perform the inversion step within inversion_type
s= inversion_type_SAR(d,A,n_patches,n_comp,inversion_opt,kkramp);
SARramp=[];
if SparseFringe_flag==1
    n_sparse_datasets=numel(SparseConstraint);
    tmp=s(end-kkramp+1:end);
    kir=1;
    for kk=1:n_sparse_datasets,
        if SparseWeight(kk)>0,
            SARramp=[SARramp;tmp(kir:kir+2)];
            kir=kir+3;
        else
            SARramp=[SARramp;0;0;0];
        end
    end
    SARramp=reshape(SARramp,3,n_sparse_datasets);
    s=s(1:end-kkramp);
end
s_var = (A'*A)^(-1);
% If positivity_flag is 1, then check if the inversion was in the correct
% orientation. This may assume a fixed direction on the time function V,
% I'm not sure.
if positivity_flag==1,
    % Currently we can only guarantee positivity if we have one monotonic time function.
    if n_comp > 1
        disp('WARNING: Positivity cannot currently be strictly imposed')
        disp('for more than one component.');
    end
    %%%%% BUILD THE STRIKE AND UPDIP SLIP VECTORS FOR EACH PATCH %%%%%
    
    sx=s(1:n_patches); sy=s(n_patches+1:2*n_patches);
    rake_curr=zeros(n_patches,1);
    rake_curr(:)=atan2(sy(:),sx(:))*180/pi;
    %%%%% CHECK IF SLIP IS IN THE SAME DIRECTION AS THE %%%%%
    %%%%% TECTONIC VECTOR TECT_VECT %%%%%
    %%%%% IF NOT, RERUN CHANGING THE SIGN OF DATA %%%%%
    
    %check positivity on all slip vectors)
    % pos_check = (abs(abs(abs(mod(rake_curr-(rake+(V_sign-1)*360),360))-180)-180)>90 || abs(sy)+abs(sx) < epsilon);
    
    % if (sum(pos_check)<numel(pos_check) || sum(abs(s)) < epsilon),
    %         disp('sense of motion NOT CONSISTENT with all tectonic vectors or inversion failed.');
    %         disp('rerunning with opposite sign ');
    %         direction_factor=-direction_factor;
    %         inversion_opt=...
    %             {positivity_flag,direction_factor,fnnls_flag,pinv_flag,...
    %             fixed_rake_flag,norm_flag,lsqnonneg_flag,rake};
    %
    %         %add lsqlin options
    %         n_inversion_opt=numel(inversion_opt);
    %         if lsqlin_options_flag,inversion_opt{n_inversion_opt+1}=lsqlin_options;end
    %
    %         s = inversion_type(d,A,n_patches,n_comp,inversion_opt);
    %         s_var = (A'*A)^(-1);
    % %        s = direction_factor*s;
    %     else
    disp('sense of motion CONSISTENT with tectonic vector');
    %     end
end
% Reshape the slip at depth vector so that G*L = U is approximately true.
L     = reshape(s,n_patches*2,n_comp);
%L_var = reshape(s_var,n_patches*2,n_comp);
L_var=s_var;
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