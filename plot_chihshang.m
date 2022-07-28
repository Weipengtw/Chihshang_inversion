%% settings
load taiwan_coastline.dat
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
p = fault_model(:,15) == 0;
ft_c = fault_model(p,1:2);
ft_c = ft_c(1:2:end,:);
xv_size = size(X_dat{1},1);
xv_e = 1:3:xv_size;
xv_n = 2:3:xv_size;
xv_u = 3:3:xv_size;

nvect = zeros(size(xv_u,2),1);
scalefactor=0.3;

%% plot slip at depth
figure;
subplot(1,2,1);
vector_scale=0.2;
test_plot_patches_point_source;
hold on

hold on

plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)


hold on
plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 80])
ylim([-40 100])
caxis([0 40])
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')

xv_h = [xv_e,xv_n];
rmsgps_h =  sqrt(sum((X_pred{1}(xv_h) - X_dat{1}(xv_h)).^2)/size(X_pred{1}(xv_h),1));
rmsgps_v = sqrt(sum((X_pred{1}(xv_u) - X_dat{1}(xv_u)).^2)/size(X_pred{1}(xv_u),1));


chi2gps_h = sum((X_pred{1}(xv_h) - X_dat{1}(xv_h)).^2./(X_err{1}(xv_h).^2))/size(X_pred{1}(xv_h),1);
chi2gps_v = sum((X_pred{1}(xv_u) - X_dat{1}(xv_u)).^2./(X_err{1}(xv_u).^2))/size(X_pred{1}(xv_u),1);

%rmsgps = rmsgps/size(X_model{1},1);
%view(27,67)
%caxis([0 70])
subplot(1,2,2);
vector_scale=0;
plot_options={'Perspective',ViewAngles,...
    'ColorMap',mycolormap,'AutoScale',1.2,'Figure','Old','ColorBar','Shading','Flat','ColorBarLabel',['SLIP (',observation_unit,')']};
plot_field_patches_point_source(fault_model,field,plot_options);view(2);hold on;
hold on
hold on
load taiwan_coastline.dat
hold on
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
%load('ps_marion_new.mat', 'lonlat_f')
%ft = lonlat_f{1};
%ft_c = llh2localxy([ft(:,2),ft(:,1)]',fliplr(origin));
hold on
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 80])
ylim([-40 100])
caxis([0 40])

RAMP = cell2mat(ramp(1));
RAMP = round(RAMP,2);
RAMP = RAMP';

set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
 title({[ 'Fixed rake  ' '\lambda_T: ',num2str(lambda_tecto),...
     '  etm0: ',num2str(etm0),]  ['  ang-tect: ',num2str(ang_tect), ...
     '  chi2: ',num2str(chi2_dense_red(1))],' V_vector: ',num2str(RAMP),},'FontWeight','normal')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


 savename=(['Joint_inversion_slip_at_depth_ang_tect_',num2str(ang_tect),'.png']);
 saveas(gcf,savename,'png');

% savename=(['test_slip_distribution_gps',num2str(X_rescale{1}),'insar',num2str(X_rescale{2}),'lambda_T',num2str(lambda_tecto),'tec_anger',num2str(ang_tect),'.png']);
% saveas(gcf,savename,'png');
%% plot horizontal velocity
figure;
subplot(1,2,1)
set(colorbar,'Visible','off');
hold on
load taiwan_coastline.dat
hold on
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
hold on
plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
axis equal

xlim([-20 80])
ylim([-40 100])


%v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_dat{1}(xv_e)).*scalefactor,(X_dat{1}(xv_n)).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_dat{1}(xv_e)-ramp{1}(1)).*scalefactor,(X_dat{1}(xv_n)-ramp{1}(2)).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');

hold on
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_model{1}(xv_e)-ramp{1}(1)).*scalefactor,(X_model{1}(xv_n)-ramp{1}(2)).*scalefactor,'LineWidth',2,'Color','r','AutoScale','off');

    
v_vector_index = quiver(60,0,40.*scalefactor,0.*scalefactor,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-5,'40 mm/yr','FontSize',10,'FontName','Hiragino Sans CNS');
v_vector_index = quiver(60,-10,40.*scalefactor,0.*scalefactor,'LineWidth',2,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-15,'data','FontSize',10,'FontName','Hiragino Sans CNS');
v_vector_index = quiver(60,-20,40.*scalefactor,0.*scalefactor,'LineWidth',2,'Color','r','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-25,'model','FontSize',10,'FontName','Hiragino Sans CNS');
vect_title = text(-15,90,'Horizontal','FontSize',25,'FontName','Hiragino Sans CNS');
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')



subplot(1,2,2) % plot residual
set(colorbar,'Visible','off');
hold on
load taiwan_coastline.dat
hold on
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
hold on
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
%plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 80])
ylim([-40 100])

%v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_model{1}(xv_e)-X_dat{1}(xv_e)).*scalefactor,(X_model{1}(xv_n)-X_dat{1}(xv_n)).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_model{1}(xv_e)-X_dat{1}(xv_e)).*scalefactor,(X_model{1}(xv_n)-X_dat{1}(xv_n)).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');

hold on
v_vector_index = quiver(60,0,40.*scalefactor,0.*scalefactor,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-5,'40 mm/yr','FontSize',10,'FontName','Hiragino Sans CNS');
vect_title = text(-15,90,'Residual','FontSize',25,'FontName','Hiragino Sans CNS');


text(-15,80,['chi2 = ',num2str(chi2gps_h)],'FontSize',15,'FontName','Hiragino Sans CNS')
text(-15,70,['rms = ',num2str(rmsgps_h) '  mm'],'FontSize',15,'FontName','Hiragino Sans CNS')

set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


savename=(['Joint_inversion_gps_horizontal_ang_tect_',num2str(ang_tect),'.png']);
saveas(gcf,savename,'png');


%savename=(['test_model_gps_horizontal_gps',num2str(X_rescale{1}),'insar',num2str(X_rescale{2}),'lambda_T',num2str(lambda_tecto),'tec_anger',num2str(ang_tect),'.png']);
%saveas(gcf,savename,'png');

%% plot Vertical
figure;
subplot(1,2,1)
set(colorbar,'Visible','off');
hold on
load taiwan_coastline.dat
hold on
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)

%load('ps_marion_new.mat', 'lonlat_f')
%ft = lonlat_f{1};

hold on
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
%plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
axis equal

xlim([-20 80])
ylim([-40 100])
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),X_dat{1}(xv_e).*0,X_dat{1}(xv_u).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');
hold on
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),X_model{1}(xv_e).*0,X_model{1}(xv_u).*scalefactor,'LineWidth',2,'Color','r','AutoScale','off');

v_vector_index = quiver(60,0,40.*scalefactor,0.*scalefactor,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-5,'40 mm/yr','FontSize',10,'FontName','Hiragino Sans CNS');
v_vector_index = quiver(60,-10,40.*scalefactor,0.*scalefactor,'LineWidth',2,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-15,'data','FontSize',10,'FontName','Hiragino Sans CNS');
v_vector_index = quiver(60,-20,40.*scalefactor,0.*scalefactor,'LineWidth',2,'Color','r','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-25,'model','FontSize',10,'FontName','Hiragino Sans CNS');
vect_title = text(-15,90,'Vertical','FontSize',25,'FontName','Hiragino Sans CNS');
%title('Horizontal')

set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')


subplot(1,2,2)
set(colorbar,'Visible','off');
hold on
load taiwan_coastline.dat
hold on
tw_c = llh2localxy([taiwan_coastline(:,2),taiwan_coastline(:,1)]',fliplr(origin));
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
load('ps_marion_new.mat', 'lonlat_f')
hold on
%plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 80])
ylim([-40 100])
v_vector = quiver(xy_position{1}(:,1),xy_position{1}(:,2),(X_model{1}(xv_e)-X_dat{1}(xv_e)).*0,(X_model{1}(xv_u)-X_dat{1}(xv_u)).*scalefactor,'LineWidth',2,'Color','k','AutoScale','off');

hold on
v_vector_index = quiver(60,0,40.*scalefactor,0.*scalefactor,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',0.8);
vect_text = text(60,-5,'40 mm/yr','FontSize',10,'FontName','Hiragino Sans CNS');
vect_title = text(-15,90,'Residual','FontSize',25,'FontName','Hiragino Sans CNS');
text(-15,80,['chi2 = ',num2str(chi2gps_v)],'FontSize',15,'FontName','Hiragino Sans CNS')
text(-15,70,['rms = ',num2str(rmsgps_v) '  mm'],'FontSize',15,'FontName','Hiragino Sans CNS')
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

savename=(['Joint_inversion_gps_vertical_ang_tect_',num2str(ang_tect),'.png']);
saveas(gcf,savename,'png');


%savename=(['test_model_gps_vertical_gps',num2str(X_rescale{1}),'insar',num2str(X_rescale{2}),'lambda_T',num2str(lambda_tecto),'tec_anger',num2str(ang_tect),'.png']);
%saveas(gcf,savename,'png');

%% plot insar

figure
subplot(1,3,1)
%test_plot_patches_point_source;
h_p = scatter(xy_position{2}(:,1),xy_position{2}(:,2),20,X_pred{2}(:),'fill');
hold on
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
hold on
plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 45])
ylim([-40 100])
colormap(jet)
caxis([-5 25])
%cb = colorbar;
%cb.Limits = [min(X_pred{1}(:))*0.9 max(X_pred{1}(:))*0.9];
set(get(colorbar,'ylabel'),'string','mm/yr');
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
vect_title = text(-15,90,'Modelled displacement','FontSize',15,'FontName','Hiragino Sans CNS');



RAMP = cell2mat(ramp(2));
RAMP = round(RAMP,2);
RAMP = RAMP';

set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
 title({[ 'Fixed rake  ' '\lambda_T: ',num2str(lambda_tecto),...
     '  etm0: ',num2str(etm0),]  ['  ang-tect: ',num2str(ang_tect), ...
     '  chi2: ',num2str(chi2_dense_red(2))],' offset: ',num2str(RAMP),},'FontWeight','normal')

subplot(1,3,2)

hold on;
[RA] = xy_position{2}(:,1) .* ramp{2}(1) + xy_position{2}(:,1) .* ramp{2}(2) +  ramp{2}(3);

h_p = scatter(xy_position{2}(:,1),xy_position{2}(:,2),20,X_dat{2}(:),'fill');
colormap(jet)
hold on
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
hold on
%plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 45])
ylim([-40 100])
colormap(jet)
caxis([-5 25])

set(get(colorbar,'ylabel'),'string','mm/yr');
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
vect_title = text(-15,90,'Observation displacement','FontSize',15,'FontName','Hiragino Sans CNS');


subplot(1,3,3)
h_p_x = scatter(xy_position{2}(:,1),xy_position{2}(:,2),20,X_pred{2}-X_dat{2},'fill');
colormap(jet)
hold on
plot(tw_c(:,1),tw_c(:,2),'Color','k','LineWidth',1)
hold on
%plot(ft_c(:,1),ft_c(:,2),'Color','k','LineWidth',1)
p = fault_model(:,9) < 1;
plot(fault_model(p,1),fault_model(p,2),'Color','k','LineWidth',1)
axis equal
xlim([-20 45])
ylim([-40 100])
%cb = colorbar;
caxis([-10 10])
colorbar
%cb.Limits = [min(X_pred{1}(:)-X_dat{1}(:))*0.9 max(X_pred{1}(:)-X_dat{1}(:))*0.9];
set(get(colorbar,'ylabel'),'string','mm/yr');
vect_title = text(-15,90,'Residual','FontSize',15,'FontName','Hiragino Sans CNS');
set(gca,'FontWeight','normal','fontsize',12,'FontName','Hiragino Sans CNS')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% Plot Ramp


savename=(['Joint_inversion_InSAR_ang_tect_',num2str(ang_tect),'.png']);
saveas(gcf,savename,'png');

% savename=(['test_model_insar_gps',num2str(X_rescale{1}),'insar',num2str(X_rescale{2}),'lambda_T',num2str(lambda_tecto),'.png']);
% saveas(gcf,savename,'png');

