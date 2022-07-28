function [G_ramp, m0,n_ramp_dataset,GF_ramp] = set_input4ramp(G, all_position, n_patches,data_type);
% Function to add a ramp for InSAR inversion and v vector for GPS
% for insar
% [d] = [ G x y 1 ] * [L a b c]' for each track
% Where: a, b and c define a ramp of equation ax+by+c
% for GPS
% [d] = [G a b c] * [L v_east v_north v_up]
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
   %G_ramp{i} = [G ramp_GF];
end
 G_ramp = []; %gps
 GF_ramp{1} = repmat(eye(3),[size(all_position{1},1),1]); %gps%gps
 GF_ramp{2} = [ramp_GF]; %insar
G_ramp = [G ramp_GF]; %insar
% Add [a b c] to the model parameters (increase size of m0 matrix)
m0 = [zeros(1,n_patches*2+3*n_dataset)];

% function [G_ramp, m0,n_ramp_dataset] = set_input4ramp(G, all_position, n_patches,data_type);
% % Function to add a ramp for InSAR inversion
% % [d] = [ G x y 1 ] * [L a b c]' for each track
% % Where: a, b and c define a ramp of equation ax+by+c
% % Warning: not tested for TS inversion
% % Written by MR
% ramp_GF = [];
% G_ramp = G;
% n_dataset = numel(all_position);
% n_ramp_dataset = numel(data_type);
% for i=1:n_dataset
%     % Add [x y 1] to the G matrix
%     ramp_GF = [zeros(length(all_position{i}),(i-1)*3) ...
%         all_position{i} ones(length(all_position{i}),1) ...
%         zeros(length(all_position{i}),(n_dataset-i-1)*3+3)];
%     % Add [a b c] to the model parameters (increase size of m0 matrix)
% end
%     G_ramp = [G_ramp ramp_GF];
% m0 = [zeros(1,n_patches*2+3*n_dataset)] ;