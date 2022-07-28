function sparse_constraint = sparse_constraint_calc(G_projected_sparse, S, V, X_time_index_sparse,n_comp)
%SPARSE_CONSTRANT_INSAR_CALC   Find the sparse constraint matrix for the
%                              inversion step
%   SPARSE_CONSTRAINT = SPARSE_CONSTRAINT_CALC(G_projected_sparse, S, V,
%   X_time_index_sparse,N_COMP) currently takes the Greens function G_projected_sparse
%   for a single sparse data set and combines it with S and V at the times from
%   X_time_index_sparse to give a linear equation for the surface displacement
%   of the sparse data set given our N_COMP model.
%
%   Example:
%   PCAIM_driver
%
%   See also PCAIM_driver.

%   By Hugo Perfettini from the code sparse_constraint_InSAR_calc.m written
%   by Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/12/08  $
%   $Revision: 1.2.0.0 $  $Date: 2010/04/27  $
%   Modified by Hugo Perfettini to treat arbitrary number of sparse sets.

n_sparse_datasets = numel(G_projected_sparse);
sparse_constraint=cell(n_sparse_datasets,1);

for kk=1:n_sparse_datasets,
    % The number of patches times the number of sources at each patch (2)
    n_patches_times_two = size(G_projected_sparse{kk},2);
    %%% Following the computations on page ??? of the manual.
    % Combine the proper entries of V with the weights expressed by S
    V_prime = S * (V(X_time_index_sparse{kk}(2),:) - V(X_time_index_sparse{kk}(1),:))';
    % Repeat this vector n_patches_times_two times for matrix computation
    V_prime_dist_prep_1 = repmat(V_prime',n_patches_times_two,1);
    % Form a diagonal matrix in preparation to multiply the columns of G_projected_sparse
    % by the proper factors.
    V_prime_dist =diag(V_prime_dist_prep_1(:));

    % construct the sparse_constraint matrix from the Greens functions and the
    % factors from V and S.
    sparse_constraint{kk} = repmat(G_projected_sparse{kk},1,n_comp) * V_prime_dist;
    
end