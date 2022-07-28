function [G_proj,G_ramp, m0, n_ramp_dataset] = set_input4ramp(G, all_position, n_patches,data_type)
% Function to add a ramp for InSAR inversion
% [d] = [ G x y 1 ] * [L a b c]' for each track
% Where: a, b and c define a ramp of equation ax+by+c
% Warning: not tested for TS inversion
% Written by MR
n_dataset = numel(all_position);
ramp_GF = [];
n_ramp_dataset = numel(data_type);
for i=1:n_dataset
    % Add [x y 1] to the G matrix
    ramp_i = [zeros(length(all_position{i}),(i-1)*3) ...
        all_position{i} ones(length(all_position{i}),1) ...
        zeros(length(all_position{i}),(n_dataset-i-1)*3+3)];
   ramp_GF = [ramp_GF; ramp_i];
   G_ramp{i} = [G ramp_GF];
end
% Add [a b c] to the model parameters (increase size of m0 matrix)
m0 = [zeros(1,n_patches*2+3*n_dataset)];
