function [X_model,X_ramp]= add_ramptomodel(X_dat,X_dat_sparse,X_model_sparse,data_type_sparse,all_position,SARramp)
%  ADD_RAMPtoMODEL
%  if the option 'Sparse_fringe' is selected in inversion parameters file the trend infered by the inversion is added
%  otherwise the programm searchs a possible trend in the residual (observed SAR data minus Modeled SAR data)
%  X_dat{i}(k,l) observation for set i, time series k, at epoch X_TIME{i}(l).
%  X_dat_sparse{i}   observation for the ieme SAR
%  X_model_sparse{i}  model for the ieme SAR
%  data_type_sparse{i} type of the data
%  all_position geographical location of the observation
%  SARramp trend infered by the inversion (ax+by+c)
%  By D. Remy
%  2010

if isempty(X_dat_sparse)
    X_model=[];
    X_ramp=[];
else
    n_dense=numel(X_dat);
    n_sparse=numel(X_model_sparse);
    sparse_data_type={'InSAR'};
    n_sparse=numel(data_type_sparse);
    for ii=1:n_sparse,
        switch data_type_sparse{ii}{1}
            case sparse_data_type,
                if isempty(SARramp)
                    disp('search a residual trend after inversion');
                    fprintf('%s \n',data_type_sparse{ii}{2});
                    long_sparse=all_position{ii+n_dense}(:,1);lat_sparse=all_position{ii+n_dense}(:,2);
                    xy_sparse=llh2localxy([lat_sparse,long_sparse]',[mean(lat_sparse),mean(long_sparse)]);
                    Xres=X_dat_sparse{ii}-X_model_sparse{ii};
                    A=[xy_sparse(:,1),xy_sparse(:,2),ones(size(xy_sparse(:,2)))];
                    ramp=A\Xres;
                    X_model{ii}=X_model_sparse{ii}+(xy_sparse(:,1)*ramp(1)+xy_sparse(:,2)*ramp(2)+ramp(3));
                    X_ramp{ii}=xy_sparse(:,1)*ramp(1)+xy_sparse(:,2)*ramp(2)+ramp(3);
                    fprintf('dx pour 100 km = %f cm dy pour 100 km = %f cm  shitf=%f \n',ramp(1)*100,ramp(2)*100,ramp(3));
                    Xres=X_dat_sparse{ii}-X_model_sparse{ii};
                    fprintf('%s residual=%8.4f instrasic rms =%8.4f\n',data_type_sparse{ii}{2},sqrt(mean(Xres.^2)),sqrt(mean(X_dat_sparse{ii}.^2)));
                    
                else
                    disp('trend infered by the  inversion');
                    fprintf('%s \n',data_type_sparse{ii}{2});
                    long_sparse=all_position{ii+n_dense}(:,1);lat_sparse=all_position{ii+n_dense}(:,2);
                    xy_sparse=llh2localxy([lat_sparse,long_sparse]',[mean(lat_sparse),mean(long_sparse)]);
                    X_model{ii}=X_model_sparse{ii}+(xy_sparse(:,1)*SARramp(1,ii)+xy_sparse(:,2)*SARramp(2,ii)+SARramp(3,ii));
                    X_ramp{ii}=xy_sparse(:,1)*SARramp(1,ii)+xy_sparse(:,2)*SARramp(2,ii)+SARramp(3,ii);
                    fprintf('dx pour 100 km = %f m dy pour 100 km = %f m  shitf=%f \n',SARramp(1,ii)*100,SARramp(2,ii)*100,SARramp(3,ii));
                    Xres=X_dat_sparse{ii}-X_model{ii};
                    fprintf('%s residual=%8.4f instrasic rms =%8.4f\n',data_type_sparse{ii}{2},sqrt(mean(Xres.^2)),sqrt(mean(X_dat_sparse{ii}.^2)));
                end
        end
    end
end
