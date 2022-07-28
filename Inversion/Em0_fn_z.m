function Em0=Em0_fn_z(fault_model,z_min,z_max,Em0_zmin,Em0_zmax,slope,options)
z=fault_model(:,3);
Em0=Em0_zmin*ones(size(z));
ind_z_sup_zmin=find(z>z_min);
% 
% 
% z_min=120;
% z_max=max(fault_model(:,3));
% dz_trans=1;
% Em0 = Em0_fn_z(fault_model,z_min,z_max,etm0,0,dz_trans,'expo');

switch options
    case 'expo' %exponential decreasing
        alpha=(Em0_zmin-Em0_zmax)/(exp(-z_min/slope)-exp(-z_max/slope));
        beta=Em0_zmin-alpha*exp(-z_min/slope);
        Em0(z<=z_min)=Em0_zmin;
        Em0(ind_z_sup_zmin)=alpha*exp(-z(ind_z_sup_zmin)/slope)+beta;
        
%     case 'Heaviside' %linear decreasing
%         Em0(ind_z_sup_zmin)=Em0_zmax;
        
    otherwise
        error(['option ',options,' is not a valid option. Choose expo or Heaviside.'])
end
