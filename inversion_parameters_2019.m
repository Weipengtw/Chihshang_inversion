
rake=fault_model(:,end);
%rake(:,1) = 100;
%p = rake < 0;
%rake = rake + 90;

n_patches=size(fault_model,1);
n_stations=size(X_dat{1},1)/3;

Aeq=[];beq=[];Ain=[];bin=[];ub=[];lb=[];s_start=[];options_lsqlin={};


%% FIXED RAKE
if numel(V)==1,
    %% VAR RAKE CASE STATIC
    signV=sign(V);
else
    %% VAR RAKE CASE DYNAMIC
    signV=sign(V(end,1)-V(1,1));
end

fault_position = fault_model(:,1:3);
area_trg=fault_model(:,6);
dl=mean(sqrt((area_trg)));

%% smoothing parameters
lambda0 = dl*ones(1,n_patches); % caracteristic patch size
%l2 = sqrt(max(area_trg)./area_trg);
%lambda0 = dl./l2';
% length scale for smoothing (km) 
lambda_tecto =40; % in the direction of tectonic rake %from where????? d
lambda_non_tecto = 40; % in the direction perpendicular tectonic rake
% if Fixed rake is used, only lambda_tecto is considered

% Standard deviation on model parameter sigma_m = etm0 = 10^(exp_etm0)
% (will later be normalized by lambda0 / lambda)
exp_etm0 = 3; % exponent (vary between -3 - 3 to propoduce L curves )
etm0=10^(exp_etm0); 
% No variation of for different patches
% Normalize sigma_m as a function of patch area
l=(area_trg);
max_1ol=(max(1./l)); % normalizing value
%Em0 = etm0 * (max(l)./l);
Em0=etm0*1./(l.*max_1ol);
%Em0(1:921) = 1000;
%Em0(922:2024) = 900;
%Em0(2025:end) = 800;
%Em0=etm0*(max(l)./l);
%Em0=etm0./(max(l)./l);
%Em0=Em0';
Em0=etm0*ones(1,n_patches);

%z_min=30;
%ind_z_sup_zmin=find(fault_model(:,12)>z_min);
%Em0(ind_z_sup_zmin) = 10^-1;
% z_max=max(fault_model(:,3));
% dz_trans=1;
%Em0 = Em0_fn_z(fault_model,z_min,z_max,etm0,10^-1,dz_trans,'expo');
%Em0 = Em0';
%Em0(ind_z_sup_zmin) = 10^-5;
which_smoothing = 'exponential'; % choices are gaussian or exponential

% find index of Y > 100km
 %idx_north = fault_model(:,3) <=10 ;
 %Em0(idx_north) = -0;
%idx_em0 = 0.01*mean(Em0);
%idx_depth = find(fault_model(:,12)>29);


% do not compute laplacian
lap_flag=0;

%% Positivity parameters
rho = 1;
n_iter =2000;

%% space selection
% 'data_space' is to be prefered, because avoid computation of
%              inv(Cm). But not implemented with positivity constraint.
% 'model_space' is better if positivity constrain is used
which_space = 'model_space'; 

%% Ramp invertion scaling
% for InSAR => ramp ax+by+c
% for cGPS => vector (a,b,c) 
% for  RES  => beta (0,0,1)
%ramp_scaling{1} = [1e-2, 1e-2, 1e5]; % 3 values for 3 ramp parameters (if we want to invert only a constant the first two values are small)
ramp_scaling{1} = [1e5 1e5 1e-5];
ramp_scaling{2} = [1e-5, 1e-5, 1e-5];
%ramp_scaling{3} = [1e-5,1e-5, 1e5];

%% Seletected options
invert_options = {...
    'Rake',rake,...
    'Vsign',signV,...
    'Uvar',U_var,...
    'WhichSmoothing',which_smoothing,...
    'WhichSpace',which_space,...
    'Lambda0',lambda0,...
    'LambdaTecto',lambda_tecto,...
    'LambdaNonTecto',lambda_non_tecto,...
    'Em0',Em0,... 
    'InvertRamp', ramp_scaling,...
    'FixedRake',...
    'Positivity',rho,n_iter,...
    };
    % 'RES_data', U_RE, G_RE,...  
   % 'InvertRamp', ramp_scaling,...
       % 'RES', m0_t...
       % 'Positivity',rho,n_iter,...
       % 'FixedRake',...
       % 'Rake',rake,...
    %'x',rho,n_iter,...

