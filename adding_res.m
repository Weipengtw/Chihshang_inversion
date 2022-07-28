
Cd = U_var;
% 
% % Set data vector
d = U(:);
G = extract_from_cell(G_projected_dense(1:2));

%% when invert best slip rate for a
%G = [G;G_projected_dense{3}];
%%
[n_tseries,n_dislocations_dirs]=size(G);
n_patches=n_dislocations_dirs/2;
n_options=size(options,2);
V_sign = 1;
%U_var = speye(numel(U));
%Em0=ones(1,n_patches);
fixed_rake_flag = 1;
L_var_flag = 0;
positivity_flag = 1;
two_faults_flag = 0;
invert_ramp = 1;
%RES_data=X_dat{3};
ramp=[];
norm_flag=1;
%Load Options




% Project the Greens function onto rake vectors
% if fixed_rake_flag ==1,  project onto rake vector
% Otherwise, projection onto 2 perpendicular directions (usefull for different smoothing in the 2 directions)

G_proj = project_Green_function(G,n_patches,fixed_rake_flag,rake,norm_flag,V_sign);

%% when invert best slip rate for a
%G_proj = [G_proj ;G_projected_dense{3}];
%%
% Set initial model concatenated

% if ramp    
[G_ramp, m0,n_ramp_dataset,GF_ramp] = set_input4ramp_peng(G_proj, all_position, n_patches, data_type);



%% add
% z=fault_model(:,12);
%  z_min = 20;
% ind_z_sup_zmin=find(z>z_min);

%m0(ind_z_sup_zmin) = sind(35)*50;
%m0(n_patches+ind_z_sup_zmin) = cosd(35)*50;
%m0(ind_z_sup_zmin) = 100;
%m0(n_patches+ind_z_sup_zmin) = 100;

%z_min=30;
%ind_z_sup_zmin=find(fault_model(:,12)>z_min);
%m0(ind_z_sup_zmin) = 40;

    
G_proj = G_ramp;
 G_proj_bk = G_proj;     

G_proj_blk = blkdiag_rep(G_proj,n_comp);
m0_blk = repmat(m0,1,n_comp);




%     m0 = zeros(1,n_patches*2);
%           
% G_proj_blk = blkdiag_rep(G_proj,n_comp);
% m0_blk = repmat(m0,1,n_comp);
% end if ramp

% compute Cm
% Set matrix of block diagonal repetition of Cm

        Cm_non_null = compute_covariance_matrix_model(fault_position,n_patches,Em0,lambda_tecto,lambda0,which_smoothing);

    Cm_null = zeros(n_patches,n_patches);
    Cm = blkdiag(Cm_non_null,Cm_null);
    Cm_blk = blkdiag_rep(Cm,n_comp);

    
    
    
    
    

%     Cm_ramp=[];
%     n_ramp = 0;
    
    
    
    % if ramp
        for i=1:n_ramp_dataset
        rs([(i-1)*3+1:(i-1)*3+3]) = ramp_scaling{i};
    end
    Cm_ramp =  diag(rs) ;
    Cm_blk = blkdiag(Cm_blk, Cm_ramp);
    n_ramp = 3*n_ramp_dataset;
    % end if ramp

% Set the initial design matrix equal to a block-diagonal repetition of the
% projected Greens function
G_proj_blk_t = G_proj_blk';



        Cm_mod = blkdiag(Cm(1:n_patches,1:n_patches), Cm_ramp);
        
        %m0_mod = zeros(1,n_patches+n_ramp);
        m0_mod = [m0(1:n_patches),zeros(1,n_ramp)];

       G_proj_mod=[G_proj(:,1:n_patches) G_proj(:,end-n_ramp+1:end); G_RE(:,1:n_patches) zeros(size(G_RE,1),n_ramp)]; % add RE Green's function
       %G_proj_mod=[G_proj(:,1:n_patches) G_proj(:,end-n_ramp+1:end)];
        

                % compute inverse model
                 m1 = admmLeastSquare_PCAIM(G_proj_mod, Cd, Cm_mod, ...
                     m0_mod, d, rho, n_iter, n_ramp);
                 disp('ADMM inversion with Positivity constrain performed')

%              Cm_blk_inv = inv(Cm_mod);
%              Gt = G_proj_mod';
%              Cd_inv = inv(Cd);  
%              A = (Gt* Cd_inv * G_proj_mod + Cm_blk_inv);
%             % compute final model
%              m1 = m0_mod' + A\Gt * Cd_inv * (d - G_proj_mod * m0_mod');
                
                
      
           m1 = [m1(1:n_patches); zeros(n_patches,1);m1(n_patches+1:end)];
        if L_var_flag
            % compute final covariance matrix % L = m0 + (Gt * Cd^{-1} * G + C_m^{-1})^{-1} * Gt * Cd^{-1} * (U - G * m0)
            Cm_1 = inv(G_proj_blk_t * Cd_inv * G_proj_blk + Cm_blk_inv);
        else
            disp('L_var is not computed. Output gives L_var=0');
        end
 
ramp = [];
% ramp
    for i = 1:n_ramp_dataset
        is = n_patches*2+(i-1)*3+1;
        ie = n_patches*2+(i-1)*3+3;
        ramp{i} =  m1(is:ie);
    end
    % Remove ramp from GF
    Cm_blk = Cm_blk(1:n_patches*2, 1:n_patches*2);
    G_proj_blk = G_proj_blk(:,1:n_patches*2);
    G_proj_blk_t = G_proj_blk';
    % split model vector between slip (m1) and ramp
    m1 = m1(1:n_patches*2); % m1 now contains only slip (no ramp)
% end ramp


%% Compute resolution matrix

G_proj_blk = [G_proj_blk;[G_RE,zeros(size(G_RE,1),size(fault_model,1))]];
G_proj_blk_t = G_proj_blk';


R = Cm_blk * G_proj_blk_t / (G_proj_blk * Cm_blk * G_proj_blk_t + Cd) * G_proj_blk;
R_diag = reshape(diag(R),n_patches*2,n_comp);
R_rest = reshape(sum(R'),n_patches*2,n_comp);
% recompose slip vector in case of fixed rake

    m1 = recompose_slip_vectors(rake,n_patches,n_comp,m1,V_sign,norm_flag);

% Reshape the slip at depth vector so that G*L = U is approximately true.
L = reshape(m1,n_patches*2,n_comp);
if L_var_flag,
    L_var = Cm_1; % CORRECT NOW
else
    L_var=0;
end

%Cm = Cm_mod;


disp('Inversion Finished');


[X_dat, final_offsets] = optimize_offsets_final(X_dat, X_err, ...
G_projected_dense, L, S, V, X_time_index);
disp('Slip model complete.');


ramp_GF=[];
 ramp{3} = [0;0;0];
 GF_ramp{3} = zeros(size(G_RE,1),3);
 for i=1:numel(G_projected_dense)
 G_projected_dense_ramp{i} = [G_projected_dense{i} GF_ramp{i}];
 L_ramp{i} = [L; ramp{i}];
 end
[X_pred,X_pred_sparse,X_model,X_model_sparse,slip_cum_ramp] = ...
create_predictions_ramp(X_dat,X_dat_sparse,G_projected_dense_ramp, ...
G_projected_sparse,L_ramp,S,V,X_time_index,X_time_index_sparse);
% Remove ramp parameters from model to prediction
slip_cum = slip_cum_ramp(1:end-numel(G_projected_dense)*3/numel(G_projected_dense));
model_statistics;

% for i=1:numel(G_projected_dense)
% L_r{i} = L;
% end
% [X_pred,X_pred_sparse,X_model,X_model_sparse,slip_cum] = ...
% create_predictions_ramp(X_dat,X_dat_sparse,G_projected_dense,...
% G_projected_sparse,L,S,V,X_time_index,X_time_index_sparse);
% model_statistics;
%plot_chihshang