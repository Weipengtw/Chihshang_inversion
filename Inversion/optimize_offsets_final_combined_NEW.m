function [X_dat,final_offsets] = optimize_offsets_final_combined_NEW(X_dat,X_err,G_projected_dense,L,S,V,X_time_index)
%OPTIMIZE_OFFSETS_FINAL   Take a dislocation model and optimize offsets
%   [X_DAT,FINAL_OFFSETS] =
%   OPTIMIZE_OFFSETS_FINAL(X_DAT,X_ERR,G,L,S,V,X_TIME_INDEX) finds the
%   constant offsets for each timeseries required to minimize final model
%   chi-squared. The amount of these offsets is FINAL_OFFSETS, and X_DAT is
%   returned with these offsets subtracted.
%
%   Example:
%   PCAIM_driver
%
%   See also invert_components, create_G_index_all, PCAIM_driver.

%   By Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/12/03  $

% load G_combined.mat;

G=extract_from_cell(G_projected_dense);

% Calculate the number of datasets
n_datasets = numel(X_time_index);

n_comp=size(S,1);

if numel(V) == 1
    for i = 1:n_datasets
        final_offsets=zeros(size(X_dat{i}));
        return;
    end
end

% Preallocate X_model, the place where we will store the model comparable
% to the original data. for comparison
X_model = cell(n_datasets,1);
for ii = 1:n_datasets
    % Calculate the number of epochs for later repetition of the means
    n_epochs_ii = size(X_dat{ii},2);
    % Calculate the model of X_dat
    %     X_model{ii} = G(X_G_index{ii},:) * L * S * V(X_time_index{ii},:)';
    
    %
    %     keyboard
    X_model_tmp=zeros(size(G*L(:,1)*V(X_time_index{ii},1)'));
    for jj=1:n_comp,
        X_model_tmp=X_model_tmp+G*L(:,jj)*S(jj,jj)*V(X_time_index{ii},jj)';
    end
    X_model{ii}=X_model_tmp;
    
    %     X_model{ii} = G * L(:,1) * S(1,1) * V(X_time_index{ii},1)'+...
    %         G * L(:,2) * S(2,2) * V(X_time_index{ii},2)';
    
    % Calculate the weighted means of both the model and data according to
    % the weights defined by X_err
    model_mean = w_mean(X_model{ii},X_err{ii},1);
    data_mean =  w_mean(X_dat{ii},X_err{ii},1);
    
    % Compute the final offsets and subtract them from X_dat
    final_offsets{ii} = data_mean - model_mean;
    X_dat{ii} = X_dat{ii} - repmat(final_offsets{ii},1,n_epochs_ii);
end