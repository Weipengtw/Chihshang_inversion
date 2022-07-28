 function [s] = recompose_slip_vectors(rake,n_patches,n_comp,s,V_sign,norm_flag)
% Recomposition of Slip Vectors in Case of a Fixed Rake (Ask APK if you
% have questions)

%if rake_flag==1,
    %%% Deal with rake of 90 mod 180
    % Find rake of 90 for a single component, and number of such patches
    rake_90_index_single_comp = find(mod(rake,180) == 90);
    rake_180_index_single_comp = find(mod(rake,180) == 0);
    n_rake_90_patches = numel(rake_90_index_single_comp);
    n_rake_180_patches = numel(rake_180_index_single_comp);
    
    % repeat the n_patches*2 offsets necessary to index the rake=90 patches
    % for all components (as each component has n_patches * 2 entries)
    repeated_slip_patch_offsets_rake_90 =repmat((0:n_comp-1)*n_patches*2,n_rake_90_patches,1);
    repeated_slip_patch_offsets_rake_180 =repmat((0:n_comp-1)*n_patches*2,n_rake_180_patches,1);
    
    % build the proper indexes for rake=90 patches' strike-slip component
    rake_90_index = repmat(rake_90_index_single_comp(:),n_comp,1) +...
        repeated_slip_patch_offsets_rake_90(:);
    rake_180_index = repmat(rake_180_index_single_comp(:),n_comp,1) +...
        repeated_slip_patch_offsets_rake_180(:);
    % set the ss component of the rake=90 patches equal to zero
    if numel(rake_90_index) > 0,
        s(rake_90_index)=0;
        
        % get the ds component of the rake=90 patches equal to the proper sign
        s(rake_90_index+n_patches) = s(rake_90_index+n_patches).*...
            repmat(sign(sind(rake(rake_90_index_single_comp(:)))),n_comp,1);
    end
    
    if numel(rake_180_index) > 0
        % set the ds component of the rake=180 patches equal to zero
        s(rake_180_index+n_patches)=0;
        % get the ss component of the rake=90 patches equal to the proper sign
        s(rake_180_index) = s(rake_180_index).*...
            repmat(sign(cosd(rake(rake_180_index_single_comp(:)))),n_comp,1);
    end
    %%% Deal with rake not equal to 90 mod 180
    % Find rake of ~90 for a single component, and number of such patches
    rake_n90_index_single_comp = find(((mod(rake,180) == 90) + (mod(rake,180) == 0))==0);
    n_rake_n90_patches = numel(rake_n90_index_single_comp);
    
    
    % repeat the n_patches*2 offsets necessary to index the rake~=90 patches
    % for all components (as each component has n_patches * 2 entries)
    repeated_slip_patch_offsets_rake_n90 =repmat((0:n_comp-1)*n_patches*2,n_rake_n90_patches,1);
    % build the proper indexes for rake~=90 patches' strike-slip component
    rake_n90_index = repmat(rake_n90_index_single_comp(:),n_comp,1) +...
        repeated_slip_patch_offsets_rake_n90(:);
    
    % set the ds component of the rake~=90 patches equal to ss component
    % times tan(rake)
    
    s(rake_n90_index+n_patches)=V_sign*s(rake_n90_index) .* ...
        repmat(sign(cosd(rake(rake_n90_index_single_comp))),n_comp,1).* ...
        tand(repmat(rake(rake_n90_index_single_comp),n_comp,1));
    s(rake_n90_index) = s(rake_n90_index).* repmat(sign(cosd(rake(rake_n90_index_single_comp))),n_comp,1);
    % If we're normalizing the slip, scale it appropriately
    if norm_flag==1,
%         disp('Normalizing Green Functions...');
        ...          s(rake_n90_index) = s(rake_n90_index) ./ repmat(sqrt(1+tand(90-rake(rake_n90_index_single_comp)).^2),n_comp,1)*sqrt(2);
            ...            s(rake_n90_index+n_patches) = s(rake_n90_index+n_patches) ./ repmat(sqrt(1+tand(rake(rake_n90_index_single_comp)).^2),n_comp,1)*sqrt(2);
            s(rake_n90_index) = s(rake_n90_index) ./ repmat(sqrt(1+tand(rake(rake_n90_index_single_comp)).^2),n_comp,1);
        s(rake_n90_index+n_patches) = s(rake_n90_index+n_patches) ./ repmat(sqrt(1+tand(rake(rake_n90_index_single_comp)).^2),n_comp,1);
    end
    
%end