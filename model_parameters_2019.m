% Greens function directory (with fortran's code)
GreensExternalFcnDir = 'Code/FaultRelated/GREENFUNC/bin/';

% Tectonic angle - to set for Fixed rake
ang_tect=145; %以正東逆時針方向旋轉
vect_tect=[cosd(ang_tect);sind(ang_tect);0];

% Fault model (computed with make_fault_model_Oaxaka.m) - Point source
% (Triangles)
% File contains 15 colums
% [xc,yc,zc,strike,dip,area,vertices] - vertices contains 9 colums
load_fault_model = [working_dir 'Chihshang_2019/fault_model/cs_37km_new4.trg'];
% fault_cs_flat_29km.trg

% Fault model - Rectangular sources
% File contains 7 columns
% [rake,area,vertices,strike_vect,updip_vect,normal_vect]
%load_fault_model = [working_dir 

% Origin file (gives lat-long of origin point)
origin_file = [working_dir 'Chihshang_2019/fault_model/origin']; 

% GF computation method
Method='Okada';

fault_file=load(load_fault_model);
strike=fault_file(:,4);strike_ang=mean(strike);
dip=fault_file(:,5);dip_angle=mean(dip);

%test_find_edges;

n_neighbours=32;

% not used
laplacian_options={'n_neighbours',n_neighbours,'free_surface_depth',0,'scaling_edge_factor',1,...
    'strike_angle',strike_ang,'dip_angle',dip_angle,'projected'};

% TRIANGULAR MODEL
get_fault_model_options = {...
    'LoadFaultModel',load_fault_model,...
    'Origin',origin_file,...
    'BuildLaplacian',laplacian_options,...
    'BuildGreensFunction',all_position, GreensExternalFcnDir,Method,...
    'TectonicVector', vect_tect};


% RECTANGULAR MODEL
% get_fault_model_options = {...
%     'LoadFaultModel',load_fault_model_file,...
%     'Origin',origin_file,...
%     'BuildLaplacian',laplacian_options,...
%     'RectangleFault',...
%     'BuildGreensFunction',all_position, GreensExternalFcnDir,...
%     'TectonicVector', vect_tect ...
% };

