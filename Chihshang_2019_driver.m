clear all;
%close all;
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');%SEEMS NECESARY FOR OSX 10.9

format long;

quite_mode_flag = 0;
addpath(genpath('Code/'));% adding main code directory
addpath(genpath('Chihshang_2019/'));% adding working directory
%working_dir='/Users/radiguem/Recherche/PCAIM_tutorial/';
working_dir='/Users/pengwei/Dropbox/kinematic_model_Chihshang/Model_Code/PCAIM/';
%working_dir='/Users/radiguem/Dropbox/kinematic_model_Chihshang/Model_Code/PCAIM/';
load_scenario_information_2019;
load_preferences;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD SCENARIO DATA FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading Data...');
run(scen_parameters_file); % open([scen_parameters_file,'.m']);


[X_time,position,X_dat,X_err,data_info,data_type]=load_all_data(data_file,first_epoch,...
    last_epoch,time_unit,sig_time,observation_unit,X_time,position,X_dat,X_err,stn_name,data_type,options);

%X_err{1}(end) = 10e3;
%X_err{1}(end-1) = 10e3;
%X_err{1}(end-2) = 10e3;

%X_err{1}(52) = 10e3;
%X_err{1}(53) = 10e3;
%X_err{1}(54) = 10e3;

%Optionnal correction of the GPS data from a constant (if we want to remove a vector and not invert it).
% diff_corr = [-15.58,23.77,0];
% % Correct GPS data from constant
%  for i=1:length(X_dat{1})/3
%      ii=(i-1)*3+1;
%      X_dat{1}(ii) = X_dat{1}(ii)-diff_corr(1);
%      X_dat{1}(ii+1) = X_dat{1}(ii+1)-diff_corr(2);
%      X_dat{1}(ii+2) = X_dat{1}(ii+2)-diff_corr(3);
%  end

 %X_dat{2} = X_dat{2} - 7;
 

X_rescale = {1,1};
sum_X_rescale=0;
for i=1:length(X_dat)
    % normalise the X_rescale so that \sum X_rescale =1
    sum_X_rescale =  sum_X_rescale + X_rescale{i};
end
for i=1:length(X_dat)
    X_rescale{i} = X_rescale{i}/sum_X_rescale;
    X_err_scaled{i} = X_err{i} .* 1/X_rescale{i};
    X_weight{i} = (1./X_err{i}).^2 * X_rescale{i};
end

% Create timeline and separate sparse data
X_time{1}=2008;
X_time{2}=2008;

[X_time_index, timeline] = create_timeline(X_time);
[X_time_index, position, X_dat, X_err, data_info, data_type, X_weight,...
    X_time_index_sparse, position_sparse, X_dat_sparse, X_err_sparse,...
    data_info_sparse,data_type_sparse, X_weight_sparse, all_position]...
    = separate_sparse_data(X_time_index, position, X_dat, X_err, ...
    data_info, data_type, X_weight, sparse_types);
disp('Data Loaded.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trick to perform static inversion.
% When static data are provided, the decomposition (Xdat = U.S.V') is not
% needed, and a single static inversion is performed. We thus set the
% U = static offsets
% S = 1
% V = 1
% A single date is considered
for i=1:length(X_time_index)
    X_time_index{i}=1;
end
% We skip the centering and decomposition and assign 1 to S (eigenvalue) V
% and V (temporal eigenvector)
% U are the static offsets
U=extract_from_cell(X_dat);
%U = [X_dat{1};X_dat{2}];
V=1;S=1;

% Compute U_var from X_err (trivial for static: diag(U_var = X_err^2))
U_var_method = 'Simple';
U_var = find_U_var(X_dat,X_err_scaled,V,S,X_time_index,n_comp,U_var_method);

%save data_InSAR_mexico_2017.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAULT MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data_InSAR_mexico_2017.mat;
disp('Building Fault model');
run(model_parameters_file);  
% open([model_parameters_file '.m'])

 load matlab_files/G.mat G;
 load matlab_files/Lap.mat Lap; 
 load matlab_files/fault_model.mat fault_model;
 load matlab_files/origin.mat origin;
 load matlab_files/rectangular_fault_flag.mat rectangular_fault_flag;
 load matlab_files/iedge.mat iedge;
 load matlab_files/G_projected_dense.mat G_projected_dense;
 load matlab_files/G_projected_sparse.mat G_projected_sparse;
 load matlab_files/sparse_constraint.mat sparse_constraint;

%   [G,Lap,fault_model,origin,rectangular_fault_flag,iedge]  = get_fault_model_NEW(get_fault_model_options);
%   [G_projected_dense,G_projected_sparse,sparse_constraint] = ...
%       project_all_greens_fcn(G,all_position,X_dat,X_dat_sparse,...
%       data_type,data_type_sparse,data_info,data_info_sparse,S,V,X_time_index_sparse,n_comp);
% % 
%        save matlab_files/G.mat G;
%        save matlab_files/Lap.mat Lap; 
%        save matlab_files/fault_model.mat fault_model;
%        save matlab_files/origin.mat origin;
%        save matlab_files/rectangular_fault_flag.mat rectangular_fault_flag;
%        save matlab_files/iedge.mat iedge;
%        save matlab_files/G_projected_dense.mat G_projected_dense;
%        save matlab_files/G_projected_sparse.mat G_projected_sparse;
%        save matlab_files/sparse_constraint.mat sparse_constraint;
% % 

for i=1:numel(all_position)
    xy_position{i} = llh2localxy([all_position{i}(:,2) all_position{i}(:,1)]',fliplr(origin));
end
disp('Fault model Complete');

%save fault_model_InSAR_mexico_2017.mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INVERT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load fault_mod\l_InSAR_mexico_2017.mat;
disp('Starting Inversion');
run(inversion_parameters_file);



    %save data_InSAR_mexico_2017.mat;
%RES=load('RES2.dat');
RES=load('RES_2022.dat');
%RES_s=load('mean_slip_repeater.txt');
RES_une=load('RES_2022_une.dat');
%load('RES_data_v2.mat','A','duration')
%load('RES_slip_rate.dat');
%load('RES_slip_rate_une.dat');

%load('RES_data_v2.mat', 'cov')

%p = cov <= 0.5;
%cov = cov(p,:);


%RES(:,6) = 0;
%RES(:,6) = RES_slip_rate(:,5);
%RES(:,7) = RES_slip_rate_une(:,5);
%RES = sortrows(RES,4);
%p = RES(:,6) < 70;
%RES = RES(p,:);
%RES(:,6) = 23;
id=3;

n_RE = size(RES,1);
G_RE = zeros(n_RE, n_patches);

for i = 1 : size(RES,1); %numbers of RES
%numbers of fault patch
d_temp = sqrt((RES(i,1) - fault_model(:,1)).^2 + (RES(i,2) - fault_model(:,2)).^2 + (RES(i,5) - fault_model(:,9)).^2);
J = find(d_temp <= 5);
J2 = find(d_temp(J)==min(d_temp(J)));
G_RE(i,J(J2)) = 1;

end

g_sum = [];
% remove G_RE with zero
for i = 1 : size(G_RE,1);
g_sum(i,1) = sum(G_RE(i,:));
end

p = find(g_sum > 0);
G_RE = G_RE(p,:);
RES= RES(p,:);
%RES_s = RES_s(p,:);
RES_une=RES_une(p,:);
%duration = duration(p,:);
%A = A(p,:);
%p = find(g_sum > 0);
%cov = cov(p);
%G_RE = G_RE./(max(sum(G_RE)));
G_projected_dense{id} = [G_RE G_RE];
X_dat{id} = RES(:,14);
%X_dat{id} = RES_s*10;
X_err{id} = (RES_une(:,8));


all_position{3} = RES(:,3:4);
position{3} = RES(:,3:4);
xy_position{3} = RES(:,1:2);
%X_rescale = {0.3,0.55,0.15};
%X_rescale = {0.4,0.6,0.01};
X_rescale = {0.4,0.45,0.15};
sum_X_rescale=0;
for i=1:length(X_dat)
% normalise the X_rescale so that \sum X_rescale =1
sum_X_rescale =  sum_X_rescale + X_rescale{i};
end
for i=1:length(X_dat)
X_rescale{i} = X_rescale{i}/sum_X_rescale;
X_err_scaled{i} = X_err{i} .* 1/X_rescale{i};
X_weight{i} = (1./X_err{i}).^2 * X_rescale{i};
end
% Create timeline and separate sparse data
X_time{1}=2008;
X_time{2}=2008;
X_time{3}=2008;
[X_time_index, timeline] = create_timeline(X_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trick to perform static inversion.
for i=1:length(X_time_index)
X_time_index{i}=1;
end
% We skip the centering and decomposition and assign 1 to S (eigenvalue) V
% and V (temporal eigenvector)
% U are the static offsets
U=extract_from_cell(X_dat);
%U = [X_dat{1};X_dat{2}];
V=1;S=1;
% Compute U_var from X_err (trivial for static: diag(U_var = X_err^2))
U_var_method = 'Simple';
U_var = find_U_var(X_dat,X_err_scaled,V,S,X_time_index,n_comp,U_var_method);

