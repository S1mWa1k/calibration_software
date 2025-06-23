close all

base_folder='D:\Cunyi_bumblebees';
date='2025_06_19';
folder='Calib_mraw';

addpath('sub_functions');

draw_images=1;

d=dir(fullfile(base_folder,date,folder,'*.cihx'));
Cameras={d.name};

num_cam=length(Cameras);

focal_length_mm=[50 50 50 50];
pixel_size=20;

fc_init = repmat(1e3.*focal_length_mm./pixel_size,2,1);
cc_init=NaN(2,num_cam);
kc_init = repmat(zeros(5,1),1,num_cam);
alpha_c=0;

dx_all=NaN(1,num_cam);
dy_all=NaN(1,num_cam);

x_cam_all=[];

% I_check=imread('new_checkerboard.tif');
% [imagePoints,boardSize] = detectCheckerboardPoints(I_check(:,:,1));
% ctrl_marker_pos=[4 4];
% save new_checkerboard_points I_check imagePoints boardSize 

load('new_checkerboard_points','boardSize','I_check','imagePoints','ctrl_marker_pos');

for c=1:num_cam

    camera=Cameras{c};

    [x_cam,XX_cam,dx,dy, nx0, ny0]=findCheckerboards(camera,base_folder,date,folder,draw_images,boardSize,I_check,imagePoints,ctrl_marker_pos);

    cc_init(:,c) = [dx;dy]/2 - 0.5; % initialize at the center of the image

    x_cam_all=[x_cam_all; x_cam];

    dx_all(c)=dx;
    dy_all(c)=dy;

    % return
end

nx=boardSize(1)-1;
ny=boardSize(2)-1;

for c=1:num_cam

    [om_grids{c}, T_grids{c}, x_cam_all(c,:), removed, residual_init_MSE, residual_init_MaxSE]=init_extrinsic_params(x_cam_all(c,:),XX_cam,fc_init(:,c),cc_init(:,c),kc_init(:,c),alpha_c);

    % any(removed)
    [fc_0dist_0cc(:,c), T_grids_0dist_0cc{c}, om_grids_0dist_0cc{c}, x_cam_all(c,:), removed, residual_0dist_0cc]=single_cam_calib_no_cc(fc_init(:,c),cc_init(:,c),kc_init(:,c),T_grids{c},om_grids{c},x_cam_all(c,:),XX_cam);
    % any(removed)

    [fc_0dist(:,c), cc_0dist(:,c), T_grids_0dist{c}, om_grids_0dist{c}, residual_0dist]=single_cam_calib_cc(fc_0dist_0cc(:,c),cc_init(:,c),kc_init(:,c),T_grids_0dist_0cc{c},om_grids_0dist_0cc{c},x_cam_all(c,:),XX_cam);

    [fc_opt(:,c), cc_opt(:,c), kc_opt(:,c), T_grids{c}, om_grids{c}, x_cam_all(c,:), removed, residual_opt]=single_cam_calib_cc_dist(fc_0dist(:,c),cc_0dist(:,c),kc_init(:,c),T_grids_0dist{c},om_grids_0dist{c},x_cam_all(c,:),XX_cam);

end

ref_cam=1;

[T_cams_2_ref, om_cams_2_ref, om_grids_2_ref, T_grids_2_ref, ref_grid, x_cam_all]=est_extrinsics_new(fc_opt,cc_opt,kc_opt,om_grids,T_grids,x_cam_all,XX_cam,'off');
[om_cams_opt,T_cams_opt,om_grids_opt,T_grids_opt,residuals]=multi_cam_calib_new(fc_opt,cc_opt,kc_opt,om_cams_2_ref,T_cams_2_ref,om_grids_2_ref,T_grids_2_ref,x_cam_all,XX_cam,num_cam,'off',nx,ny,ref_grid,ref_cam);

RR=hypot(residuals(1:2:end,:),residuals(2:2:end,:));
MSE=mean(RR(:),'omitnan');

disp(['Mean error: ' num2str(MSE) 'px'])

calib_data.fc=fc_opt;
calib_data.cc=cc_opt;
calib_data.kc=kc_opt;
calib_data.x_cam_all=x_cam;
calib_data.XX_cam=XX_cam{1};
calib_data.num_cam=num_cam;
calib_data.om_cams=om_cams_opt;
calib_data.T_cams=T_cams_opt;
calib_data.om_grids=om_grids_opt;
calib_data.T_grids=T_grids_opt;
calib_data.residuals=residuals;
calib_data.ref_cam=ref_cam;
calib_data.cam_names=Cameras;
calib_data.nx=nx;
calib_data.ny=ny;

[cam_distance, cam_angle]=plot_calibration_setup_function(calib_data);

save(fullfile(base_folder,date,folder,'calibration_file'),'calib_data')
