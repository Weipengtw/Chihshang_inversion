function sparse_constraint = sparse_constraint_InSAR_calc(G_InSAR, S, V, InSAR_time_index,n_comp)
%SPARSE_CONSTRANT_INSAR_CALC   Find the sparse constraint matrix for the
%                              inversion step
%   SPARSE_CONSTRAINT = SPARSE_CONSTRANT_INSAR_CALC(G_InSAG_INSAR, S, V,
%   INSAR_TIME_INDEX,N_COMP) currently takes the Greens function G_INSAR
%   for a single InSAR image and combines it with S and V at the times from
%   INSAR_TIME_INDEX to give a linear equation for the surface displacement
%   of the InSAR image given our N_COMP model.
%
%   Example:
%   PCAIM_driver
%
%   See also PCAIM_driver.

%   By Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/12/08  $
%   $Revision: 1.2.0.0 $  $Date: 2010/04/27  $
%   Modified by Hugo Perfettini to treat arbitrary number of sparse sets.

n_sparse_datasets = numel(G_InSAR);
sparse_constraint=cell(n_sparse_datasets,1);

for kk=1:n_sparse_datasets,
    % The number of patches times the number of sources at each patch (2)
    n_patches_times_two = size(G_InSAR{kk},2);
    %%% Following the computations on page ??? of the manual.
    % Combine the proper entries of V with the weights expressed by S
    if numel(V)>1,
        V_prime = S * (V(InSAR_time_index{kk}(2),:) - V(InSAR_time_index{kk}(1),:))';
    else
        V_prime=S*V;%When there is a unique measure (interseismic, coseismic,..)
    end
    % Repeat this vector n_patches_times_two times for matrix computation
    V_prime_dist_prep_1 = repmat(V_prime',n_patches_times_two,1);
    % Form a diagonal matrix in preparation to multiply the columns of G_InSAR
    % by the proper factors.
    V_prime_dist =diag(V_prime_dist_prep_1(:));
    
    % construct the sparse_constraint matrix from the Greens functions and the
    % factors from V and S.
    sparse_constraint{kk} = repmat(G_InSAR{kk},1,n_comp) * V_prime_dist;
end