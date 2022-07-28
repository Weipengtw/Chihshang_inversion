function s=inversion_type(d,A,n_patches,n_comp,inversion_opt)
%INVERSION_TYPE   Perform the actual inversion operation given options
%   S=inversion_type(D,A,N_PATCHES,N_COMP,INVERSION_OPT) gives the best
%   solution to the inverse problem A*S = D given the number of components
%   N_COMP and the options in INVERSION_OPT.
%
%   Example:
%   PCAIM_driver
%
%   See also INVERT_COMPONENTS, PCAIM_DRIVER.
%

%   Written and commented by Hugo Perfettini
%   Additional comments and revised by Andrew Kositsky
%   Rewritten to reduce run time
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.1.0.0 $  $Date: 2010/06/30  $

% Set the options, assuming this order
positivity_flag=inversion_opt{1};
fact=inversion_opt{2};
fnnls_flag=inversion_opt{3};
pinv_flag=inversion_opt{4};
rake_flag=inversion_opt{5};
norm_flag=inversion_opt{6};
lsqnonneg_flag=inversion_opt{7};
rake=inversion_opt{8};
V_sign = inversion_opt{9};

Aeq=[];beq=[];Ain=[];bin=[];s_start=[];options_lsqlin={};
lb=[zeros(2*n_patches,1);-inf*ones((n_comp-1)*2*n_patches,1)];
ub=zeros(2*n_patches + (n_comp-1)*2*n_patches,1)+Inf;
lsqlin_options={};
linprog_options={};
lsqlin_flag=0;linprog_flag=0;

if numel(inversion_opt)>9
    [lsqlin_flag_tmp,ind_lsqlin]=find_string_in_cell(inversion_opt{10},'lsqlin');
    if lsqlin_flag_tmp || isempty(inversion_opt{10}), lsqlin_options=inversion_opt{10};lsqlin_flag=1;end
    [linprog_flag,ind_linprog]=find_string_in_cell(inversion_opt{10},'linprog');
    if linprog_flag, linprog_options=inversion_opt{10};end
end

if linprog_flag==1 && lsqlin_flag==1, error('lsqlin and linprog can not be called simultaneously');end

%loading lsqlin options if needed
if lsqlin_flag,
    [iflag,index]=find_string_in_cell(lsqlin_options,'lb');if iflag,lb=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'ub');if iflag,ub=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'Aeq');if iflag,Aeq=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'beq');if iflag,beq=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'Ain');if iflag,Ain=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'bin');if iflag,bin=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'s_start');if iflag,s_start=lsqlin_options{index+1};end
    [iflag,index]=find_string_in_cell(lsqlin_options,'options_lsqlin');if iflag,options_lsqlin=lsqlin_options{index+1};end
end

if linprog_flag,
    [iflag,index]=find_string_in_cell(linprog_options,'f');if iflag,f=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'lb');if iflag,lb=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'ub');if iflag,ub=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'Aeq');if iflag,Aeq=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'beq');if iflag,beq=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'Ain');if iflag,Ain=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'bin');if iflag,bin=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'s_start');if iflag,s_start=linprog_options{index+1};end
    [iflag,index]=find_string_in_cell(linprog_options,'options_linprog');if iflag,options_linprog=linprog_options{index+1};end
    linprog_flag=1;
end

%Inversion with Positivity Constrain

if positivity_flag==1,
    if fnnls_flag==1 && lsqnonneg_flag==0,
        disp('check positivity...');
        if lsqlin_flag==1,
            disp('lsqlin used');
            % Note that when there are more parameters we're solving for
            % (s) than in lb or ub, lb and ub are applied only to the first
            % numel(lb) parameters of s.
            %             s=lsqlin(A,fact*d,Ain,bin,Aeq,beq,lb,ub,s_start,options_lsqlin);
            [s,resnorm,residual,exitflag,output,lambda]=lsqlin(A,fact*d,Ain,bin,Aeq,beq,lb,ub,s_start,options_lsqlin);
            
            message=['exitflag: ',num2str(exitflag)];disp(message);
%             if exitflag~=1, error('!!! Problem with lsqlin. Consider changing the parameters (e.g., smoothing)');end
            
        elseif linprog_flag==1,
            % Note that when there are more parameters we're solving for
            % (s) than in lb or ub, lb and ub are applied only to the first
            % numel(lb) parameters of s.
            %             s=linprog(f,Ain,bin,A,fact*d,lb,ub,s_start,options_linprog);
            [s,fval,exitflag,output,lambda] = linprog(f,Ain,bin,A,fact*d,lb,ub,s_start,options_linprog);            
            keyboard
            message=['exitflag: ',num2str(exitflag)];disp(message);
            %              if exitflag~=1, error('!!! Problem with linprog. Consider changing the parameters (e.g., smoothing)');end
        else
            disp('WARNING: lsqlin not present (most likely, the optimization toolbox');
            disp('         has not been installed).');
            XtX=A'*A;Xty=A'*d;
            disp('fnnls used');
            %
            % fnnls assumes that ALL components must have positive slip.
            % This is generally NOT the case for multicomponent models.
            %
            s=fnnls(XtX,fact*Xty);
        end
    end
    if lsqnonneg_flag==1,
        disp('lsqnonneg used');
        XtX=A'*A;Xty=A'*d;
        s=lsqnonneg(A,fact*d);
    end
else
    %Inversion without Positivity Constrain
    if pinv_flag==1,
        % Inversion using Pseudo Inverse
        disp('pseudo inverse used');
        s=pinv(A) * d;
    else
        % Classical Least-Squares Inversion
        s=A\d;
    end
end

% plot(d,'r.');
% hold on;
% plot(A*s,'bo');
% keyboard

% Recomposition of Slip Vectors in Case of a Fixed Rake (Ask APK if you
% have questions)
if rake_flag==1,
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
    
end

% Recomposition of Slip Vectors in case of Positivity Constrain
if positivity_flag==1,
    s=fact.*s;
end