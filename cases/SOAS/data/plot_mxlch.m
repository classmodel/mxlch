function plot_mxlch(mxlch_path)
% plot_mxlch is used to load the observations and model outputs and
% perform visualization on the dataset.
% Inputs:
%     mxlch_path: path of the MXLCH output folder
close all; % close all figures
tzos = -6; % time zone offset: CST = UTC - 6 h
save_figure = false; % save figure or not

%--------------------------------------------------------------------------
% Load diurnal VOCs from top of AABC tower
% f_voc = 'D:\SOAS_2013\MXLCH_data\Chem\VOC\AABC_tower_averaged.mat';
% vocs_aabc_diurnal = load(f_voc,'voc_avg');
% vocs_aabc_diurnal = vocs_aabc_diurnal.voc_avg;
% vocs_aabc_diurnal_avg = vocs_aabc_diurnal.voc_avg_ppbv;
% vocs_aabc_diurnal_std = vocs_aabc_diurnal.voc_std_ppbv;
% vocs_aabc_diurnal_spt = vocs_aabc_diurnal.sampling_time_middle_hours;
% % Find the index range for each chunk of data separated by NaN
% ind_data = zeros(20,1);
% pt_ind = 1; % pointer for the ind_data
% pt_this = 1; % pointer for voc data
% ind_data(pt_ind) = pt_this;
% pt_ind = pt_ind + 1;
% voc_len = size(vocs_aabc_diurnal_avg,1);
% for i=2:voc_len
%     if isnan(vocs_aabc_diurnal_avg(i)) ~= isnan(vocs_aabc_diurnal_avg(pt_this))
%         ind_data(pt_ind) = i - 1;
%         pt_this = i;
%         pt_ind = pt_ind + 1;
%         ind_data(pt_ind) = pt_this;
%         pt_ind = pt_ind + 1;
%     end
% end
% ind_data(pt_ind) = voc_len;
% ind_data = ind_data(1:pt_ind);
% ind_data = reshape(ind_data,2,[]);
% ind_data_diff = diff(ind_data);
% ind_data(:,ind_data_diff==0) = [];
% num_block = size(ind_data,2);
%--------------------------------------------------------------------------
% Load O3, NOx, VOCs from C130, only use the data from June 12, 2013
f_C130 = 'D:\SOAS\Master Thesis\C130_merged.mat';
d_C130 = load(f_C130);
d_C130 = d_C130.rf04; % only use data from June 12, 2013
% Set the altitude range to 100-1000 m
ind_use = d_C130.Altitude > 100 & d_C130.Altitude < 1000;

var_n = {'theta_K', 'SH_g_per_kg', 'O3_ppbv', 'NO_pptv', 'NO2_pptv', ...
         'ISOP_ppbv', 'MVK_ppbv', 'MT_ppbv'};
for i = 1 : numel(var_n)
    C130.([var_n{i}, '_avg']) = nanmean(d_C130.(var_n{i})(ind_use));
    C130.([var_n{i}, '_std']) = nanstd(d_C130.(var_n{i})(ind_use));
end

rf04_spt_UTC = d_C130.spt_UTC(ind_use);
rf04_spt_CST     = rf04_spt_UTC(round((1 + end)/2)) - 6/24;
rf04_spt_CST_hrs = datevec(rf04_spt_CST);
rf04_spt_CST_hrs = rf04_spt_CST_hrs(4:6) * [1;1/60;1/3600];

C130_spt_CST_hrs = rf04_spt_CST_hrs;
%--------------------------------------------------------------------------
% Load VOCs flux data from top of AABC tower
f_VOCs_flux = 'D:\SOAS\Master Thesis\Flux_avg.mat';
flux_data = load(f_VOCs_flux);
%--------------------------------------------------------------------------
% Load averaged VOCs from WASP system and the VOCs from AABC tower during
% the same time period
f_VOCs_WASP = 'D:\SOAS\Master Thesis\WASP_AABC_aligned.mat';
vocs_wasp = load(f_VOCs_WASP,'wasp');
vocs_wasp = vocs_wasp.wasp;
vocs_wasp_voc_avg = vocs_wasp.voc_avg;
vocs_wasp_voc_std = vocs_wasp.voc_std;
vocs_wasp_hrs = datevec(vocs_wasp.time_middle_UTC);
vocs_wasp_days = vocs_wasp_hrs(:,3); % days
vocs_wasp_hrs = tzos + vocs_wasp_hrs(:,4) + vocs_wasp_hrs(:,5)/60 + vocs_wasp_hrs(:,6)/3600;
% Ignore the cloudy days
vocs_wasp_ignore = vocs_wasp_days == 5 | vocs_wasp_days == 6 | ...
    vocs_wasp_days == 8 | vocs_wasp_days == 9; % cloudy days will be ignored

vocs_wasp_voc_avg = vocs_wasp_voc_avg(~vocs_wasp_ignore,:);
vocs_wasp_voc_std = vocs_wasp_voc_std(~vocs_wasp_ignore,:);
vocs_wasp_hrs     = vocs_wasp_hrs(~vocs_wasp_ignore,:);

% VOCs data from 06/01 RF1 and RF2
vocs_wasp_06_01_spt_CST  = [5.3;6.95];
vocs_wasp_06_01_isop     = [0.6;0.7];
vocs_wasp_06_01_isop_std = [0.19;0.22];
vocs_wasp_06_01_m71      = [0.64;0.46];
vocs_wasp_06_01_m71_std  = [0.17;0.13];
vocs_wasp_06_01_mt       = [1.04;1.12];
vocs_wasp_06_01_mt_std   = [0.19;0.18];
% Data from AABC tower which correspond to WASP sampling period
% vocs_aabc_single_point = load(f_VOCs_WASP,'aabc');
% vocs_aabc_single_point = vocs_aabc_single_point.aabc;
% vocs_aabc_single_point_avg = vocs_aabc_single_point.voc_avg;
% vocs_aabc_single_point_std = vocs_aabc_single_point.voc_std;
% vocs_aabc_single_point_spt = vocs_aabc_single_point.voc_sampling_time_UTC;
% spt_temp = datevec(vocs_aabc_single_point_spt(:,2));
% vocs_aabc_single_point_spt = tzos + spt_temp(:,4) + spt_temp(:,5)./60 + spt_temp(:,6)./3600;
% % Ignore cloudy days
% vocs_aabc_single_point_avg = vocs_aabc_single_point_avg(~vocs_wasp_ignore,:);
% vocs_aabc_single_point_std = vocs_aabc_single_point_std(~vocs_wasp_ignore,:);
% vocs_aabc_single_point_spt = vocs_aabc_single_point_spt(~vocs_wasp_ignore,:);
%--------------------------------------------------------------------------
% Load the sensible and latent heat flux data
f_met = 'D:\SOAS\Master Thesis\AABC_Sonic_flux_for_MXLCH.mat';
d_met = load(f_met,'output_all');
d_met = d_met.output_all;
%--------------------------------------------------------------------------
% Load boundary layer height data from CTR ceilometer
f_blh_ceilometer = 'D:\SOAS\Master Thesis\SEARCH_Ceilometer.mat';
d_blh_ceilometer = load(f_blh_ceilometer,'output');
d_blh_ceilometer = d_blh_ceilometer.output;
temp      = datevec(d_blh_ceilometer.sampling_time_CST);
d_blh_spt = temp(:,4) + temp(:,5)./60 + temp(:,6)./3600;
d_blh_avg = d_blh_ceilometer.BLH_avg_m;
d_blh_std = d_blh_ceilometer.BLH_std_m;
%--------------------------------------------------------------------------
% Load isoprene hydroxynitrates
f_ISOPN = 'D:\SOAS\Master Thesis\SEARCH_ISOPN.mat';
d_ISOPN = load(f_ISOPN, 'd_out');
d_ISOPN = d_ISOPN.d_out;
%--------------------------------------------------------------------------
% Load boundary layer height data from sounding
% f_blh_souding = 'D:\SOAS_2013\MXLCH_data\Boundary_layer_height\iss3_sounding\output.mat';
% d_blh_souding = load(f_blh_souding);
%--------------------------------------------------------------------------
% Load NOx, O3, HOx, and HCHO mixing ratio data from SEARCH site
f_NOx_O3 = 'D:\SOAS\Master Thesis\CTR_Gases_130715_final.mat';
NOx_O3   = load(f_NOx_O3); % O3, NO, NO2
NOx_O3   = NOx_O3.output_averaged;

gas_tt   = NOx_O3.hours; % time
O3_mean  = NOx_O3.O3_ppb;
O3_std   = NOx_O3.O3_ppb_std;
NO_mean  = NOx_O3.NO_ppb;
NO_std   = NOx_O3.NO_ppb_std;
NO2_mean = NOx_O3.NO2_ppb;
NO2_std  = NOx_O3.NO2_ppb_std;
HNO3_mean = NOx_O3.HNO3_ppb;
HNO3_std = NOx_O3.HNO3_ppb_std;
PAN_mean = NOx_O3.PAN_ppb;
PAN_std = NOx_O3.PAN_ppb_std;
%--------------------------------------------------------------------------
f_HOx  = 'D:\SOAS\Master Thesis\HOx_SEARCH_site.mat';
HOx    = load(f_HOx,'output'); % load HOx data
HOx    = HOx.output;
%--------------------------------------------------------------------------
f_HCHO = 'D:\SOAS\Master Thesis\HCHO_SEARCH_site.mat';
HCHO   = load(f_HCHO); % load HCHO data
HCHO   = HCHO.output;
%--------------------------------------------------------------------------
% Load model output
if exist(fullfile(mxlch_path,'chem_conc'),'file')==2 % chemistry is simulated
    %[dyn, sca, chem, flux] = load_MXLCH(mxlch_path);
    [dyn, sca, chem] = load_MXLCH(mxlch_path);
    h_chem = chem{1}; % header with chemical names
    d_chem = chem{2};
    d_chem(:,1) = d_chem(:,1) + tzos; % time zone offset
    %d_flux = flux{2};
    %d_flux(:,1) = d_flux(:,1) + tzos; % time zone offset
else % chemistry is not simulated
    [dyn, sca] = load_MXLCH(mxlch_path);
    d_chem = [];
end
d_dyn      = dyn{2};
d_dyn(:,1) = d_dyn(:,1) + tzos;
d_sca      = sca{2};
d_sca(:,1) = d_sca(:,1) + tzos;

% C130_t = [10.5, 10.5];
t = d_met.H_avg(:,1); % time zone offset
t = datevec(t);
t = t(:,4) + t(:,5)/60 + t(:,6)/3600; % convert to hours number
%--------------------------------------------------------------------------
rho_Cp  = 1.231e3;
rho_Lv  = 3013.5;

    
% Load flux tower average mixing ratios
% file_flux = 'D:\SOAS_2013\Flux\Flux_diurnal_avg.mat';
% flux_avg = load(file_flux);
% flux_avg = flux_avg.output;
%
% file_flux_flux = 'D:\SOAS_2013\Sonic\EC_calulation\Flux_avg.mat';
% flux_data = load(file_flux_flux);
%--------------------------------------------------------------------------
% START TO PLOT
shades_1 = [204,204,255]./[255,255,255]; % shade color for the 1 sd
shades_2 = [204,255,204]./[255,255,255]; % shade color for the 1 sd
color_wasp = 'b'; % color of wasp data

h_ind = 1; % figure handle index
% plot w't'
h(h_ind) = figure('Color','w');
% errorbar(t,d.H_avg(2:end,2),d.H_std(2:end,2),'bo');
d_low = d_met.H_avg(:,2) - d_met.H_std(:,2);
d_high = 2*d_met.H_std(:,2);
flag = ~(isnan(d_low) & isnan(d_high));
h_a = area(t(flag),[d_low(flag),d_high(flag)],-0.1,'linestyle','none');
set(h_a(1),'FaceColor',[1,1,1]);
set(h_a(2),'FaceColor',shades_1);
set(h_a,'LineWidth',0.01)

hold on;
plt_h1 = plot(t,d_met.H_avg(:,2),'b-','linewidth',2);
ylim([-0.05, 0.2])
 yyaxis right;
   ylim([-0.05, 0.2] * rho_Cp);
   ylabel('Sensible heat flux [W m^{-2}]');
   yyaxis left;
hold on;
plt_h2 = plot(d_dyn(:,1),d_dyn(:,8),'r-','linewidth',2);
grid on;
xlabel('Time (CST)');
ylabel('Sensible heat flux (K m s^{-1})');
xlim([6,18])
ylim([-0.05, 0.2])
set(gca,'XTick',0:2:24)
legend([plt_h1,plt_h2],{'AABC tower','MXLCH'});
set(gca,'layer','top');
set(findall(gcf,'type','text'),'fontSize',12)
set(gca,'FontSize',12)
h_ind = h_ind + 1;


% plot w'q'
h(h_ind) = figure('Color','w');
% errorbar(t,d.E_avg(2:end,2),d.E_std(2:end,2),'bo');
d_low = d_met.E_avg(:,2) - d_met.E_std(:,2);
d_high = 2*d_met.E_std(:,2);
flag = ~(isnan(d_low) & isnan(d_high));
h_a = area(t(flag),[d_low(flag),d_high(flag)],-0.2,'linestyle','none');

set(h_a(1),'FaceColor',[1,1,1]);
set(h_a(2),'FaceColor',shades_1);
set(h_a,'LineWidth',0.01)

hold on;
plt_h1 = plot(t,d_met.E_avg(:,2),'b-','linewidth',2);
ylim([-0.1, 0.3])
 yyaxis right;
   ylim([-0.1, 0.3] * rho_Lv);
   ylabel('Latent heat flux [W m^{-2}]');
   yyaxis left;
hold on;
plt_h2 = plot(d_sca(:,1),d_sca(:,7),'r-','linewidth',2);
grid on;
xlabel('Time (CST)');
ylabel('Latent heat flux (g kg^{-1} m s^{-1})');
xlim([6,18])
ylim([-0.1, 0.3])
set(gca,'XTick',0:2:24)
legend([plt_h1,plt_h2],{'AABC tower','MXLCH'});
set(gca,'layer','top');
set(findall(gcf,'type','text'),'fontSize',12)
set(gca,'FontSize',12)
h_ind = h_ind + 1;

% plot theta
h(h_ind) = figure('Color','w');
% errorbar(t,d.theta_avg(2:end,2),d.theta_std(2:end,2),'bo');
h_a = area(t,[d_met.theta_avg(:,2) - d_met.theta_std(:,2),2*d_met.theta_std(:,2)],290,'linestyle','none');
set(h_a(1),'FaceColor',[1,1,1]);
set(h_a(2),'FaceColor',shades_1);
set(h_a,'LineWidth',0.01)

hold on;
plt_h1 = plot(t,d_met.theta_avg(:,2),'b-','linewidth',2);

hold on;
plt_h2 = plot(d_dyn(:,1),d_dyn(:,5),'r-','linewidth',2);
grid on;
hold on;
% data from C130
h_e1 = errorbar(C130_spt_CST_hrs,C130.theta_K_avg,C130.theta_K_std);
set(h_e1                            , ...
    'LineStyle'       , 'none'      , ...
    'Color'           , [.3 .3 .3]  , ...
    'LineWidth'       , 2           , ...
    'Marker'          , 'o'         , ...
    'MarkerSize'      , 10           , ...
    'MarkerEdgeColor' , [.2 .2 .2]  , ...
    'MarkerFaceColor' , [.7 .7 .7]  );

xlabel('Sampling time (CST)');
ylabel('\theta (K)');
xlim([0,24])
set(gca,'XTick',0:2:24)
lh = legend([plt_h1,plt_h2,h_e1],{'AABC tower','MXLCH','C-130 June 12'},'location','NorthWest');
M  = findobj(lh,'type','line');
set(M,'linewidth',2) % this will also work on a vector of handles
set(gca,'layer','top');
set(findall(gcf,'type','text'),'fontSize',12)
set(gca,'FontSize',12)
h_ind = h_ind + 1;

% plot q
h(h_ind) = figure('Color','w');
% errorbar(t,d.q_avg(2:end,2),d.q_std(2:end,2),'bo');

d_low = d_met.q_avg(:,2) - d_met.q_std(:,2);
d_high = 2*d_met.q_std(:,2);
flag = ~(isnan(d_low) & isnan(d_high));
h_a = area(t(flag),[d_low(flag),d_high(flag)],12,'linestyle','none');

set(h_a(1),'FaceColor',[1,1,1]);
set(h_a(2),'FaceColor',shades_1);
set(h_a,'LineWidth',0.01)

hold on;
plt_h1 = plot(t,d_met.q_avg(:,2),'b-','linewidth',2);
hold on;
plt_h2 = plot(d_sca(:,1),d_sca(:,4),'r-','linewidth',2);
grid on;

% data from C130
h_e1 = errorbar(C130_spt_CST_hrs,C130.SH_g_per_kg_avg,C130.SH_g_per_kg_std);

set(h_e1                            , ...
    'LineStyle'       , 'none'      , ...
    'Color'           , [.3 .3 .3]  , ...
    'LineWidth'       , 2           , ...
    'Marker'          , 'o'         , ...
    'MarkerSize'      , 10           , ...
    'MarkerEdgeColor' , [.2 .2 .2]  , ...
    'MarkerFaceColor' , [.7 .7 .7]  );


xlabel('Sampling time (CST)');
ylabel('Speficif humidity (g kg^{-1})');
xlim([0,24])
set(gca,'XTick',0:2:24)
lh = legend([plt_h1,plt_h2,h_e1],{'AABC tower','MXLCH','C-130 June 12'},'location','SouthWest');
M  = findobj(lh,'type','line');
set(M,'linewidth',2) % this will also work on a vector of handles
set(gca,'layer','top');
set(findall(gcf,'type','text'),'fontSize',12)
set(gca,'FontSize',12)
h_ind = h_ind + 1;

% plot boundary layer height
% sounding data of 4,11,12,13
ind_blh = 1; % 06/04
tt{ind_blh} = [9;15] + [0;11]./60;
bdl{ind_blh} = [332;1538];

ind_blh = ind_blh + 1; % 06/11
tt{ind_blh} = [9;15] + [0;0]./60;
bdl{ind_blh} = [515;1574];

ind_blh = ind_blh + 1; % 06/12
tt{ind_blh} = [9;14] + [0;55]./60;
bdl{ind_blh} = [1119;1435];

ind_blh = ind_blh + 1; % 06/13
tt{ind_blh} = [9;15] + [0;0]./60;
bdl{ind_blh} = [537;1512];

% ind_blh = ind_blh + 1; % 06/14
% tt{ind_blh} = [9;14] + [0;54]./60;
% bdl{ind_blh} = [782;1390];
% 
% ind_blh = ind_blh + 1; % 06/15
% tt{ind_blh} = [9;14] + [0;56]./60;
% bdl{ind_blh} = [613;1803];
% 
% ind_blh = ind_blh + 1; % 06/25
% tt{ind_blh} = [9;15] + [10;1]./60;
% bdl{ind_blh} = [899;1334];

bdl = cell2mat(bdl);
bdl_avg = mean(bdl,2);
bdl_std = std(bdl,0,2);
% clrs = [255,0,0;0,255,0;0,0,255;255,0,255;128,0,0;128,128,0;0,128,0;128,0,128]./repmat([255,255,255],8,1);
h(h_ind) = figure('Color','w');

d_blh_high = 2*d_blh_std;
d_blh_low = d_blh_avg - d_blh_std;
flag = ~(isnan(d_blh_low) & isnan(d_blh_high));
axh_1 = area(d_blh_spt(flag),[d_blh_low(flag),d_blh_high(flag)],0,'linestyle','none');

set(axh_1(1),'FaceColor',[1,1,1]);
set(axh_1(2),'FaceColor',shades_1);
set(axh_1,'LineWidth',0.01)

hold on;
axh_2 = plot(d_blh_spt,d_blh_avg,'b-','linewidth',2);
hold on;
axh_3 = errorbar(tt{1},bdl_avg,bdl_std);
set(axh_3                            , ...
    'LineStyle'       , 'none'      , ...
    'Color'           , 'g'  , ...
    'LineWidth'       , 2           , ...
    'Marker'          , 's'         , ...
    'MarkerSize'      , 10           , ...
    'MarkerEdgeColor' , 'g'  , ...
    'MarkerFaceColor' , 'none'  );
hold on;
axh_4 = plot(d_dyn(:,1),d_dyn(:,3),'r-','linewidth',2);
grid on;
xlabel('Sampling time (CST)');
ylabel('Boundary layer height (m)');
xlim([6,18])
set(gca,'XTick')
set(gca,'layer','top');
% lgd_str = {'MXLCH','06/04','06/11','06/12','06/13','06/14','06/15','06/25','06/26'};
lgd_str = {'Ceilometer','Sounding','MXLCH'};
lh = legend([axh_2,axh_3,axh_4],lgd_str,'location','NorthWest');
M  = findobj(lh,'type','line');
set(M,'linewidth',2) % this will also work on a vector of handles
set(findall(gcf,'type','text'),'fontSize',12)
set(gca,'FontSize',12)
h_ind = h_ind + 1;


% plot chemicals
if ~isempty(d_chem)
    % Ozone
    O3_ind = strcmp(h_chem,'O3');
    h(h_ind) = figure('Color','w');
    d_low = O3_mean - O3_std;
    d_high = 2*O3_std;
    flag = ~(isnan(d_low) & isnan(d_high));
    h_a = area(gas_tt(flag),[d_low(flag),d_high(flag)],-10,'linestyle','none');
    
    set(h_a(1),'FaceColor',[1,1,1]);
    set(h_a(2),'FaceColor',shades_1);
    set(h_a,'LineWidth',0.01)
    
    hold on;
    h_O3_meas = plot(gas_tt,O3_mean,'b-','linewidth',2);
    hold on;
    h_O3_mxlch = plot(d_chem(:,1),d_chem(:,O3_ind),'r-','linewidth',2);
    hold on;
    
    h_e1 = errorbar(C130_spt_CST_hrs,C130.O3_ppbv_avg,C130.O3_ppbv_std);

    set(h_e1                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , [.3 .3 .3]  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , [.2 .2 .2]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    lh = legend([h_O3_meas,h_O3_mxlch,h_e1],{'SEARCH tower','MXLCH','C-130 June 12'},'location','SouthEast');
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('O3 (ppbv)');
    xlim([6,18])
    set(gca,'XTick',0:2:24)
    ylim([0,70]);
    grid on;
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % NO
    h(h_ind) = figure('Color','w');
    NO_ind = strcmp(h_chem,'NO');
    d_low = NO_mean - NO_std;
    d_high = 2*NO_std;
    flag = ~(isnan(d_low) & isnan(d_high));
    h_a = area(gas_tt(flag),[d_low(flag),d_high(flag)],-0.5,'linestyle','none');
    
    set(h_a(1),'FaceColor',[1,1,1]);
    set(h_a(2),'FaceColor',shades_1);
    set(h_a,'LineWidth',0.01)
    
    hold on;
    h_NO_meas = plot(gas_tt,NO_mean,'b-','linewidth',2);
    hold on;
    h_NO_mxlch = plot(d_chem(:,1),d_chem(:,NO_ind),'r-','linewidth',2);
    
    h_e1 = errorbar(C130_spt_CST_hrs,C130.NO_pptv_avg/1000,C130.NO_pptv_std/1000);

    set(h_e1                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , [.3 .3 .3]  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , [.2 .2 .2]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    lh = legend([h_NO_meas,h_NO_mxlch,h_e1],{'SEARCH tower','MXLCH','C130 June 12'});
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('NO (ppbv)');
    xlim([6,18])
    set(gca,'XTick',0:2:24)
    ylim([0,0.6]);
    grid on;
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % Plot NO flux data
%     h(h_ind) = figure('Color','w');
%     plot(d_flux(:,1),d_flux(:,NO_ind),'r-','linewidth',2);
%     grid on;
%     xlabel('Sampling time (CST)');
%     ylabel('NO flux (ppbv m s^{-1})');
%     xlim([0,24])
%     set(gca,'XTick',0:2:24)
%     h_ind = h_ind + 1;
    
    % NO2
    % add up all the AN listed in the new isoprene mechanism
%     AN_spc = {'ISOPNB', 'ISOPND', 'MACRN', 'MVKN', ...
%         'PAN', 'PMN', 'PPN', 'R4N2', 'PROPNN', 'ETHLN', 'MNO3'};
%     AN_mr = zeros(size(d_chem, 1), numel(AN_spc));
%     AN_ind = false(1, numel(h_chem));
%     for ind_this = 1 : numel(AN_spc)
%         AN_ind(:) = strcmp(h_chem, AN_spc{ind_this});
%         AN_mr(:, ind_this) = d_chem(:, AN_ind);
%     end
    
    h(h_ind) = figure('Color','w');
    NO2_ind = strcmp(h_chem,'NO2');
    d_low = NO2_mean - NO2_std;
    d_high = 2*NO2_std;
    flag = ~(isnan(d_low) & isnan(d_high));
    h_a = area(gas_tt(flag),[d_low(flag),d_high(flag)],-1,'linestyle','none');
    
    set(h_a(1),'FaceColor',[1,1,1]);
    set(h_a(2),'FaceColor',shades_1);
    set(h_a,'LineWidth',0.01)
    
    hold on;
    h_NO2_meas = plot(gas_tt,NO2_mean, 'b-', 'linewidth', 2);
    hold on;
    h_NO2_mxlch = plot(d_chem(:,1), d_chem(:,NO2_ind) ,'r-','linewidth',2);
    
    h_e1 = errorbar(C130_spt_CST_hrs(1),C130.NO2_pptv_avg/1000,C130.NO2_pptv_std/1000);

    set(h_e1                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , [.3 .3 .3]  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , [.2 .2 .2]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    lh = legend([h_NO2_meas,h_NO2_mxlch,h_e1],{'SEARCH tower','MXLCH','C-130 June 12'});
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('NO_{2} (ppbv)');
    xlim([6,18])
    ylim([0,2.5])
    set(gca,'XTick',0:2:24)
    grid on;
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % Plot AN
%     h(h_ind) = figure('Color','w');
%     
%     clr_use = distinguishable_colors(size(AN_mr, 2));
%     for ind_this = 1 : size(AN_mr, 2)
%         plot(d_chem(:,1), AN_mr(:, ind_this), 'linewidth', 2, 'color', clr_use(ind_this, :));
%         hold on
%     end
%     hold on
%     plot(d_chem(:,1), sum(AN_mr, 2), 'b--', 'linewidth', 2);
%     
%     legend(AN_spc);
%     xlabel('Sampling time (CST)');
%     ylabel('AN (ppbv)');
%     grid on;
%     set(findall(gcf,'type','text'),'fontSize',12)
%     set(gca,'FontSize',12)
%     h_ind = h_ind + 1;
    
    % Plot PAN
%     h(h_ind) = figure('Color','w');
%     PAN_ind = strcmp(h_chem,'PAN');
%     d_low = PAN_mean - PAN_std;
%     d_high = 2*PAN_std;
%     flag = ~(isnan(d_low) & isnan(d_high));
%     h_a = area(gas_tt(flag),[d_low(flag),d_high(flag)],-0.5,'linestyle','none');
%     
%     set(h_a(1),'FaceColor',[1,1,1]);
%     set(h_a(2),'FaceColor',shades_1);
%     set(h_a,'LineWidth',0.01)
%     
%     hold on;
%     h_PAN_meas = plot(gas_tt,PAN_mean,'b-','linewidth',2);
%     hold on;
%     h_PAN_mxlch = plot(d_chem(:,1),d_chem(:,PAN_ind),'r-','linewidth',2);
%     
% 
%     lh = legend([h_PAN_meas,h_PAN_mxlch],{'SEARCH tower','MXLCH'});
%     M  = findobj(lh,'type','line');
%     set(M,'linewidth',2) % this will also work on a vector of handles
%     xlabel('Sampling time (CST)');
%     ylabel('PAN (ppbv)');
%     xlim([0,24])
%     set(gca,'XTick',0:2:24)
%     %ylim([-0.4,1.4]);
%     grid on;
%     set(gca,'layer','top');
%     set(findall(gcf,'type','text'),'fontSize',12)
%     set(gca,'FontSize',12)
%     h_ind = h_ind + 1;
    
    
    % Plot CH2O
    CH2O_ind = strcmp(h_chem,'CH2O');
    h(h_ind) = figure('Color','w');
    plot(HCHO.hhmm_mid,HCHO.HCHO_avg,'b-','linewidth',2);
    hold on;
    plot(d_chem(:,1),d_chem(:,CH2O_ind),'r-','linewidth',2);
    legend({'SEARCH tower','MXLCH'});
    xlabel('Sampling time (CST)');
    ylabel('HCHO (ppbv)');
    grid on;
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % Plot isoprene
    iso_ind = strcmp(h_chem,'ISO');
    h(h_ind) = figure('Color','w');
    %     % AABC tower measurements
    %     d_low = vocs_aabc_diurnal_avg(:,4) - vocs_aabc_diurnal_std(:,4);
    %     d_high = 2.*vocs_aabc_diurnal_std(:,4);
    %     for ind_temp=1:num_block
    %         ind_use = ind_data(1,ind_temp):ind_data(2,ind_temp);
    %         h_a = area(vocs_aabc_diurnal_spt(ind_use),[d_low(ind_use),d_high(ind_use)],-2,'linestyle','none');
    %         set(h_a(1),'FaceColor',[1,1,1]);
    %         set(h_a(2),'FaceColor',shades_1);
    %         hold on;
    %     end
    %     h_aabc = plot(vocs_aabc_diurnal_spt,vocs_aabc_diurnal_avg(:,4),'b-','linewidth',2);
    %     hold on;
    
    h_mxlch = plot(d_chem(:,1),d_chem(:,iso_ind),'r-','linewidth',2);
    hold on;
    h_e1 = errorbar([vocs_wasp_hrs;vocs_wasp_06_01_spt_CST],[vocs_wasp_voc_avg(:,4);vocs_wasp_06_01_isop],[vocs_wasp_voc_std(:,4);vocs_wasp_06_01_isop_std]);
    set(h_e1                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , color_wasp  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , color_wasp  , ...
        'MarkerFaceColor' , 'b'  );
    hold on;
    % C130 data
    h_e2 = errorbar(C130_spt_CST_hrs,C130.ISOP_ppbv_avg,C130.ISOP_ppbv_std);

    set(h_e2                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , [.3 .3 .3]  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , '^'         , ...
        'MarkerSize'      , 10          , ...
        'MarkerEdgeColor' , [.2 .2 .2]  , ...
        'MarkerFaceColor' , [.2 .2 .2]  );

    %     h_e2 = errorbar(vocs_aabc_single_point_spt,vocs_aabc_single_point_avg(:,4),vocs_aabc_single_point_std(:,4));
    %     set(h_e2                            , ...
    %         'LineStyle'       , 'none'      , ...
    %         'Color'           , [.3 .3 .3]  , ...
    %         'LineWidth'       , 2           , ...
    %         'Marker'          , 'o'         , ...
    %         'MarkerSize'      , 8           , ...
    %         'MarkerEdgeColor' , 'r'  , ...
    %         'MarkerFaceColor' , [.7 .7 .7]  );
    %     legend([h_aabc,h_mxlch,h_e1,h_e2,h_e3],{'AABC tower','MXLCH','WASP','C130 June 12','C130 June 14'},'location','northwest');
    lh = legend([h_mxlch,h_e1,h_e2],{'MXLCH','WASP','C130 June 12'},'location','southeast');
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Time (CST)');
    ylabel('Isoprene [ppb]');
    xlim([5,18])
    set(gca,'XTick',0:2:24)
    grid on;
    xlim([5,18]);
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % Plot MVK
    mvk_ind = strcmp(h_chem,'MVK');
    macr_ind = strcmp(h_chem,'MACR');
    if any(macr_ind) % for complex mechanism, MVK and MACR are separated
        mvk_macr = d_chem(:,mvk_ind) + d_chem(:,macr_ind);
    else
        mvk_macr = d_chem(:,mvk_ind);
    end
    h(h_ind) = figure('Color','w');
    %     % AABC measurements
    %     d_low = vocs_aabc_diurnal_avg(:,5) - vocs_aabc_diurnal_std(:,5);
    %     d_high = 2.*vocs_aabc_diurnal_std(:,5);
    %     for ind_temp=1:num_block
    %         ind_use = ind_data(1,ind_temp):ind_data(2,ind_temp);
    %         h_a = area(vocs_aabc_diurnal_spt(ind_use),[d_low(ind_use),d_high(ind_use)],0,'linestyle','none');
    %         set(h_a(1),'FaceColor',[1,1,1]);
    %         set(h_a(2),'FaceColor',shades_1);
    %         hold on;
    %     end
    %     h_aabc = plot(vocs_aabc_diurnal_spt,vocs_aabc_diurnal_avg(:,5),'b-','linewidth',2); % AABC diurnal
    %     hold on;
    
    h_mxlch = plot(d_chem(:,1),mvk_macr,'r-','linewidth',2);
    hold on;
    h_e1 = errorbar([vocs_wasp_hrs;vocs_wasp_06_01_spt_CST],[vocs_wasp_voc_avg(:,5);vocs_wasp_06_01_m71],[vocs_wasp_voc_std(:,5);vocs_wasp_06_01_m71_std]);
    set(h_e1                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , color_wasp  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , color_wasp  , ...
        'MarkerFaceColor' , shades_1  );
    hold on;
    % C130 data
    h_e2 = errorbar(C130_spt_CST_hrs,C130.MVK_ppbv_avg,C130.MVK_ppbv_std);

    set(h_e2                            , ...
        'LineStyle'       , 'none'      , ...
        'Color'           , [.3 .3 .3]  , ...
        'LineWidth'       , 2           , ...
        'Marker'          , '^'         , ...
        'MarkerSize'      , 10           , ...
        'MarkerEdgeColor' , [.2 .2 .2]  , ...
        'MarkerFaceColor' , [.7 .7 .7]  );

    %     h_e2 = errorbar(vocs_aabc_single_point_spt,vocs_aabc_single_point_avg(:,5),vocs_aabc_single_point_std(:,5));
    %     set(h_e2                            , ...
    %         'LineStyle'       , 'none'      , ...
    %         'Color'           , [.3 .3 .3]  , ...
    %         'LineWidth'       , 2           , ...
    %         'Marker'          , 'o'         , ...
    %         'MarkerSize'      , 8           , ...
    %         'MarkerEdgeColor' , 'r'  , ...
    %         'MarkerFaceColor' , [.7 .7 .7]  );
    %     legend([h_aabc,h_mxlch,h_e1,h_e2,h_e3],{'AABC tower','MXLCH','WASP','C130 June 12','C130 June 14'},'location','northwest');
    lh = legend([h_mxlch,h_e1,h_e2],{'MXLCH','WASP','C130 June 12'},'location','southeast');
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('MVK+MACR (ppbv)');
    xlim([0,24])
    set(gca,'XTick',0:2:24)
    grid on;
    xlim([0,24]);
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
%     % Plot MT
%     apin_ind = strcmp(h_chem,'APIN');
%     Bpin_ind = strcmp(h_chem,'BPIN');
%     limo_ind = strcmp(h_chem,'LIMO');
%     if any(limo_ind) % for complex mechanism, MVK and MACR are separated
%        mt = d_chem(:,apin_ind) + d_chem(:,Bpin_ind) + d_chem(:,limo_ind);
%     else
%         mt = d_chem(:,apin_ind);
%     end
    
    h(h_ind) = figure('Color','w');
    %     % AABC measurements
    %     d_low = vocs_aabc_diurnal_avg(:,7) - vocs_aabc_diurnal_std(:,7);
    %     d_high = 2.*vocs_aabc_diurnal_std(:,7);
    %     for ind_temp=1:num_block
    %         ind_use = ind_data(1,ind_temp):ind_data(2,ind_temp);
    %         h_a = area(vocs_aabc_diurnal_spt(ind_use),[d_low(ind_use),d_high(ind_use)],-0.1,'linestyle','none');
    %         set(h_a(1),'FaceColor',[1,1,1]);
    %         set(h_a(2),'FaceColor',shades_1);
    %         hold on;
    %     end
    %
    %     h_aabc = plot(vocs_aabc_diurnal_spt,vocs_aabc_diurnal_avg(:,7),'b-','linewidth',2);
    %     hold on;
    
%     h_mxlch = plot(d_chem(:,1),d_chem(:,mt),'r-','linewidth',2);
%     hold on;
%     h_e1 = errorbar([vocs_wasp_hrs;vocs_wasp_06_01_spt_CST],[vocs_wasp_voc_avg(:,7);vocs_wasp_06_01_mt],[vocs_wasp_voc_std(:,7);vocs_wasp_06_01_mt_std]);
%     set(h_e1                            , ...
%         'LineStyle'       , 'none'      , ...
%         'Color'           , color_wasp  , ...
%         'LineWidth'       , 2           , ...
%         'Marker'          , 'o'         , ...
%         'MarkerSize'      , 10           , ...
%         'MarkerEdgeColor' , color_wasp  , ...
%         'MarkerFaceColor' , shades_1  );
%     hold on;
%     % C130 data
%     h_e2 = errorbar(C130_spt_CST_hrs,C130.MT_ppbv_avg,C130.MT_ppbv_std);
%     set(h_e2                            , ...
%         'LineStyle'       , 'none'      , ...
%         'Color'           , [.3 .3 .3]  , ...
%         'LineWidth'       , 2           , ...
%         'Marker'          , '^'         , ...
%         'MarkerSize'      , 10           , ...
%         'MarkerEdgeColor' , [.2 .2 .2]  , ...
%         'MarkerFaceColor' , [.7 .7 .7]  );

    %     h_e2 = errorbar(vocs_aabc_single_point_spt,vocs_aabc_single_point_avg(:,7),vocs_aabc_single_point_std(:,7));
    %     set(h_e2                            , ...
    %         'LineStyle'       , 'none'      , ...
    %         'Color'           , [.3 .3 .3]  , ...
    %         'LineWidth'       , 2           , ...
    %         'Marker'          , 'o'         , ...
    %         'MarkerSize'      , 8           , ...
    %         'MarkerEdgeColor' , 'r'  , ...
    %         'MarkerFaceColor' , [.7 .7 .7]  );
    %     legend([h_aabc,h_mxlch,h_e1,h_e2,h_e3],{'AABC tower','MXLCH','WASP','C130 June 12','C130 June 14'},'location','north');
    
%     lh = legend([h_mxlch,h_e1,h_e2],{'MXLCH','WASP','C130 June 12'},'location','NorthEast');
%     M  = findobj(lh,'type','line');
%     set(M,'linewidth',2) % this will also work on a vector of handles
%     xlabel('Time (CST)');
%     ylabel('Monoterpenes (ppb)');
%     xlim([6,18])
%     set(gca,'XTick',0:2:24)
%     grid on;
%     xlim([6,18]);
%     ylim([0,1.5]);
%     set(gca,'layer','top');
%     set(findall(gcf,'type','text'),'fontSize',12)
%     set(gca,'FontSize',12)
%     h_ind = h_ind + 1;
%     
    %     % Plot CO
    %     CO_ind = strcmp(h_chem,'CO');
    %     h(h_ind) = figure('Color','w');
    %     plot(gas_tt,CO_mean,'b-','linewidth',2);
    %     hold on;
    %     plot(d_chem(:,1),d_chem(:,CO_ind),'r-','linewidth',2);
    %     legend({'SEARCH tower','MXLCH'},'location','NorthEast');
    %     xlabel('Sampling time (CST)');
    %     ylabel('CO (ppbv)');
    %     grid on;
    %     h_ind = h_ind + 1;
%     
%     
    %Plot the flux of isoprene
%     h(h_ind) = figure('Color','w');
%     h_e1 = errorbar(flux_data.hhmm_all,flux_data.w_VOC_avg(:,4),flux_data.w_VOC_std(:,4));
%     set(h_e1                            , ...
%        'LineStyle'       , 'none'      , ...
%        'Color'           , [.3 .3 .3]  , ...
%        'LineWidth'       , 2           , ...
%        'Marker'          , 'o'         , ...
%        'MarkerSize'      , 8           , ...
%        'MarkerEdgeColor' , [.2 .2 .2]  , ...
%        'MarkerFaceColor' , [.7 .7 .7]  );
%     %remove the errorbar end ticks
%     hh_e = findall(h_e1);
%     dataLen = length( get(h_e1, 'xdata') );
%     x_org = get(hh_e(3),'xdata');
%     y_org = get(hh_e(3),'ydata');
%     ind_use = logical(kron(ones(1,dataLen),[ones(1,3),zeros(1,6)]));
%     x_new = x_org(ind_use);
%     y_new = y_org(ind_use);
%     set(hh_e(3),'xdata',x_new)
%     set(hh_e(3),'ydata',y_new)
%     
%    hold on;
%     plot(d_flux(:,1),d_flux(:,iso_ind),'r-','linewidth',2);
%     grid on;
%     lh = legend({'AABC tower','MXLCH'});
%     M  = findobj(lh,'type','line');
%     set(M,'linewidth',2) % this will also work on a vector of handles
%     xlabel('Time CST (h)');
%     ylabel('Isoprene (ppbv m s^{-1})');
%     xlim([7,18])
%     set(gca,'XTick')
%     xlim([7,18])
%     set(gca,'XTick')
%     set(findall(gcf,'type','text'),'fontSize',12)
%     set(gca,'FontSize',12)
%     h_ind = h_ind + 1;
%     
    % Plot the flux of mt
   % h(h_ind) = figure('Color','w');
   % h_e1 = errorbar(flux_data.hhmm_all,flux_data.w_VOC_avg(:,9),flux_data.w_VOC_std(:,9));
   % set(h_e1                            , ...
   %     'LineStyle'       , 'none'      , ...
       % 'Color'           , [.3 .3 .3]  , ...
       % 'LineWidth'       , 2           , ...
       % 'Marker'          , 'o'         , ...
       % 'MarkerSize'      , 8           , ...
       % 'MarkerEdgeColor' , [.2 .2 .2]  , ...
       % 'MarkerFaceColor' , [.7 .7 .7]  );
    %hold on;
    
    % remove the errorbar end ticks
   % hh_e = get(h_e1,'children');
   % dataLen = length( get(h_e1, 'xdata') );
   % x_org = get(hh_e(2),'xdata');
    %y_org = get(hh_e(2),'ydata');
    %ind_use = logical(kron(ones(1,dataLen),[ones(1,3),zeros(1,6)]));
    %x_new = x_org(ind_use);
    %y_new = y_org(ind_use);
    %set(hh_e(2),'xdata',x_new)
    %set(hh_e(2),'ydata',y_new)
    
   % plot(d_flux(:,1),d_flux(:,mt_ind),'r-','linewidth',2);
   % grid on;
   % lh = legend({'AABC tower','MXLCH'});
   % M  = findobj(lh,'type','line');
   % set(M,'linewidth',2) % this will also work on a vector of handles
   % xlabel('Time CST (h)');
   % ylabel('Monoterpenes (ppbv m s^{-1})');
   % xlim([0,24])
   % set(gca,'XTick',0:2:24)
   % xlim([0,24])
   % set(gca,'XTick',0:2:24)
   % set(findall(gcf,'type','text'),'fontSize',12)
   % set(gca,'FontSize',12)
   % h_ind = h_ind + 1;
    
    
    % Plot OH
    h(h_ind) = figure('Color','w');

    % Convert the unit from ppbv to molecules cm-3
    P = 100510; % pressure, unit = Pa
    % Load the dynamics file to get the potential temperature
    spt_theta = datevec(d_met.theta_avg(:,1));
    % convert to hours number, used for interpolation
    spt_theta = spt_theta(:,4) + spt_theta(:,5)./60 + spt_theta(:,5)./3600;
    theta = d_met.theta_avg(:,2); % obtain the potential temperature
    % Calculate the number density of air
    n_air = 6.02e23 .* P ./ (8.314 .* theta); % unit=molecules m-3
    n_air = n_air .* 1e-6; % unit=molecules cm-3
    % interpolate
    spt_HOx = HOx.sampling_time_CST; % already in hours number
    n_air_1 = interp1(spt_theta,n_air,spt_HOx);
    HOx.HOx_avg(:,1) = HOx.HOx_avg(:,1) .* 1e-9 .* n_air_1(:); % convert from ppbv to molecules cm-3
    HOx.HOx_std(:,1) = HOx.HOx_std(:,1) .* 1e-9 .* n_air_1(:); % convert from ppbv to molecules cm-3
    
    d_low = HOx.HOx_avg(:,1) - HOx.HOx_std(:,1);
    d_high = 2*HOx.HOx_std(:,1);
    flag = ~(isnan(d_low) & isnan(d_high));
    h_a = area(HOx.sampling_time_CST(flag),[d_low(flag),d_high(flag)],-0.1,'linestyle','none');
    set(h_a(1),'FaceColor',[1,1,1]);
    set(h_a(2),'FaceColor',shades_1);
    set(h_a,'LineWidth',0.01)
    
    hold on;
    h_plt1 = plot(HOx.sampling_time_CST,HOx.HOx_avg(:,1),'b-','linewidth',2);
    hold on;
    
    oh_ind = strcmp(h_chem,'OH');
    n_air_2 = interp1(spt_theta,n_air,d_chem(:,1));
    d_chem(:,oh_ind) = d_chem(:,oh_ind) .* 1e-9 .* n_air_2(:);
    
    h_plt2 = plot(d_chem(:,1),d_chem(:,oh_ind),'r-','linewidth',2);
    grid on;
    lh = legend([h_plt1,h_plt2],{'SEARCH tower','MXLCH'});
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('OH (molec. cm^{-3})');
    xlim([6,18])
    set(gca,'XTick',0:2:24)
    set(gca,'XTick',0:2:24)
    %ylim([0,18e-5]);
    set(gca,'layer','top');
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    % Plot HO2
    h(h_ind) = figure('Color','w');
    
    HOx.HOx_avg(:,3) = HOx.HOx_avg(:,3) .* 1e-9 .* n_air_1(:); % convert from ppbv to molecules cm-3
    HOx.HOx_std(:,3) = HOx.HOx_std(:,3) .* 1e-9 .* n_air_1(:); % convert from ppbv to molecules cm-3
    
    h_e1 = plot(HOx.sampling_time_CST,HOx.HOx_avg(:,3),'bo','linewidth',2);
    hold on;

    HO2_ind = strcmp(h_chem,'HO2');
    
    HO2_mxlch = d_chem(:,HO2_ind) .* 1e-9 .* n_air_2(:);
  
    plot(d_chem(:,1),HO2_mxlch,'r-','linewidth',2);
    grid on;
    lh = legend({'SEARCH tower','MXLCH'});
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('HO_2 (molec. cm^{-3})');
    xlim([0,24])
    set(gca,'XTick',0:2:24)
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;
    
    NO_HO2 = d_chem(:, NO_ind)./d_chem(:, HO2_ind); % NO:HO2

    % Plot ISOPOOH 
    h(h_ind) = figure('Color','w');

    ISOPOOH_ind = strcmp(h_chem, 'ISOPOO');
    ISOPOOH = d_chem(:,ISOPOOH_ind);
    
    plot(d_chem(:,1), ISOPOOH,'r-','linewidth',2);
    grid on;
    lh = legend('MXLCH');
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('ISOPOOH (ppbv)');
    xlim([0,24])
    set(gca,'XTick',0:2:24)
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    h_ind = h_ind + 1;

    % Plot isoprene hydroxynitrates
    h(h_ind) = figure('Color','w');

    ISOPNB_ind = strcmp(h_chem,'ISOPNB');
    ISOPND_ind = strcmp(h_chem,'ISOPND');
    
    d_low = d_ISOPN.ISOPN_avg_pptv - d_ISOPN.ISOPN_std_pptv;
    d_high = 2 * d_ISOPN.ISOPN_std_pptv;
    h_a = area(d_ISOPN.Sampling_time_CST, [d_low, d_high], -100, 'linestyle', 'none');
    set(h_a(1),'FaceColor',[1,1,1]);
    set(h_a(2),'FaceColor',shades_1);
    set(h_a,'LineWidth',0.01)
    c
    hold on;
    h_ISOPN_obs = plot(d_ISOPN.Sampling_time_CST, d_ISOPN.ISOPN_avg_pptv, 'b-','linewidth',2);
    
    ISOPN = (d_chem(:,ISOPNB_ind) + d_chem(:,ISOPND_ind)) * 1000;
    
    h_ISOPN_mxlch = plot(d_chem(:,1), ISOPN, 'r-', 'linewidth', 2);
    grid on;
    lh = legend([h_ISOPN_obs, h_ISOPN_mxlch], {'SEARCH', 'MXLCH'});
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('ISOPN (pptv)');
    xlim([0,24])
%     ylim([0, 150])
    set(gca,'XTick',0:2:24)
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    set(gca,'layer','top');
    h_ind = h_ind + 1;
    
    % Plot ISOPOOH/(MVK+MACR)
    h(h_ind) = figure('Color','w');
    plot(d_chem(:,1), ISOPOOH./mvk_macr, 'b-', 'linewidth', 2);
    grid on;
    lh = legend('MXLCH');
    M  = findobj(lh,'type','line');
    set(M,'linewidth',2) % this will also work on a vector of handles
    xlabel('Sampling time (CST)');
    ylabel('ISOPOOH/(MVK+MACR)');
    xlim([0,24])
    set(gca,'XTick',0:2:24)
    set(findall(gcf,'type','text'),'fontSize',12)
    set(gca,'FontSize',12)
    set(gca,'layer','top');
    h_ind = h_ind + 1; 
end
% Export the figures
if save_figure
    cd('D:\SOAS_2013\MXLCH_data\Figure_export');
    isOpen  = exportToPPTX();
    if ~isempty(isOpen),
        % If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    exportToPPTX('new'); % export figures to pptx
    for i=1:length(h)
        %         set(h(i), 'PaperPosition', [0 0 8 5]);
        %         set(h(i), 'PaperSize', [8 5]); %Keep the same paper size
        %saveas(h(i), ['fig_',num2str(i),'.png'], 'png');
        exportToPPTX('addslide');
        exportToPPTX('addpicture',h(i));
    end
    save_time = datestr(now,'yyyy_mm_dd_HH_MM');
    exportToPPTX('saveandclose',['MXLCH_output_',save_time]); % save and close file
end

end
%==========================================================================
function varargout = load_MXLCH(datapath)
% Load the output of MXLCH model run

% Load the file manually
f_dyn = 'output_dyn'; % dynamics output
f_sca = 'output_sca'; % humidity output
f_chem = 'chem_conc'; % chemistry conc
f_flux = 'chem_flux';

if nargout==2
    varargout{1} = get_data(fullfile(datapath, f_dyn));
    varargout{2} = get_data(fullfile(datapath, f_sca));
elseif nargout==3
    varargout{1} = get_data(fullfile(datapath, f_dyn));
    varargout{2} = get_data(fullfile(datapath, f_sca));
    varargout{3} = get_data(fullfile(datapath, f_chem));
elseif nargout==4
    varargout{1} = get_data(fullfile(datapath, f_dyn));
    varargout{2} = get_data(fullfile(datapath, f_sca));
    varargout{3} = get_data(fullfile(datapath, f_chem));
    varargout{4} = get_data(fullfile(datapath, f_flux));
else
    error('nargout should be 2 or 3!')
end

end
%==========================================================================
function output = get_data(file_name)
%% Obtain data from file
fid = fopen(file_name);
h1  = fgetl(fid); % read first header, title of this data
h2  = fgetl(fid); % read second header, number of rows
h3  = strtrim(fgetl(fid)); % read third header, header of data WITHOUT THE FUCKING WHITESPACE YOU MOTHERFUCKER

num_row  = str2double(h2); % total number of rows
header   = strsplit(h3); % split the headers
num_itms = numel(header); % how many columns
fmt      = repmat('%f',1,num_itms); % construct format for fscanf

d        = fscanf(fid, fmt, num_row*num_itms);
fclose(fid);
n        = size(d,1); % number of total data elements
if mod(n,num_itms)~=0
    error(['Number of elements mismatch! (' num2str(n) ' vs ' num2str(num_itms) ')']);
end
d        = reshape(d,num_itms,n/num_itms)'; % reshape d
output   = {header, d};


end






