function [om_cams_opt,T_cams_opt,om_grids_opt,T_grids_opt,residuals]=multi_cam_calib_new(fc,cc,kc,om_cams_2_ref,T_cams_2_ref,om_grids_2_ref,T_grids_2_ref,x_cam,XX_cam,num_cam,disp,nx,ny,ref_grid,ref_cam)
% % % 
% save temp_multi_cam_calib_0dist
% return
% clear
% close all

% load temp_multi_cam_calib_0dist
% load temp_multi_cam_calib_new

% fc=repmat(mean(fc),2,1);
% 
% return

% for each grid point, find the number of camera views its been found in
% and remove any with no camera views (a single camera view should now be
% ok).

num_grids=size(om_grids_2_ref,2);

for i=1:num_grids
   
    xx=cell2mat(x_cam(:,i));
    X=XX_cam{i};
    
    xnan=sum(~isnan(xx));
    
    xx(:,xnan<4)=[];
    X(:,xnan<4)=[];
    
    x_cell{i}=xx;
    X_cell{i}=X;
    
    grid_l(i)=size(xx,2);
    
end

% set reference grid to 0 as it won't be used in calibration
grid_l(ref_grid)=0;

% remove any grids where no points have been found

x_cell(:,grid_l==0)=[];
X_cell(grid_l==0)=[];

T_grids_2_ref(:,grid_l==0)=[];
om_grids_2_ref(:,grid_l==0)=[];
grid_l(grid_l==0)=[];

% % reduce numbers, just for testing
% use_grids=1:10;
% grid_l=grid_l(use_grids);
% 
% T_grids_2_ref=T_grids_2_ref(:,use_grids);
% om_grids_2_ref=om_grids_2_ref(:,use_grids);
% 
% x_cell=x_cell(use_grids);
% X_cell=X_cell(use_grids);
% 
% for i=1:length(x_cell)
%     x_cell{i}=x_cell{i}(:,1:5:40);
%     X_cell{i}=X_cell{i}(:,1:5:40);
%     grid_l(i)=size(x_cell{i},2);
%  
% end

% return
xx=cell2mat(x_cell);

chi0=[om_cams_2_ref T_cams_2_ref om_grids_2_ref T_grids_2_ref];

grid_t=numel(xx);

for i=1:num_cam
    col_start=3*i-2;
    
    jcol_cam_temp=repmat([col_start:col_start+2 3*(num_cam)+col_start:3*(num_cam)+col_start+2],sum(grid_l),1);
    jcol_cam_temp(isnan(xx(2*i,:)),:)=[];
    
    jcol_cam{i,1}=repmat(jcol_cam_temp,2,1);

    jrow_cam_temp=(2*i-1:num_cam*2:sum(grid_t))';
    jrow_cam_temp(isnan(xx(2*i,:)))=[];
    
    jrow_cam{i,1}=repmat([jrow_cam_temp;jrow_cam_temp+1],1,6);
    
end

jcol_cam=cell2mat(jcol_cam);
jrow_cam=cell2mat(jrow_cam);

jcol_grid=[];
jrow_grid=[];
jrow_max=0;

for i=1:length(grid_l)
    
    jcol_grid=[jcol_grid;repmat([6*(num_cam)+3*i-2:6*(num_cam)+3*i 6*(num_cam)+3*i-2+3*length(grid_l):6*(num_cam)+3*i+3*length(grid_l)],num_cam*2*grid_l(i),1)];
    jrow_grid=[jrow_grid;repmat((1:num_cam*2*grid_l(i))',1,6)+jrow_max];
    jrow_max=max(jrow_grid(:));
end

a=ismember(jrow_grid,jrow_cam);
jcol_grid(a==0)=[];
jrow_grid(a==0)=[];

jcol=[jcol_cam(:);jcol_grid(:)];
jrow=[jrow_cam(:);jrow_grid(:)];

JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,6*(num_cam)+3*length(grid_l)+3*length(grid_l));

% disp='off';

options=optimset('lsqnonlin');
options=optimset(options,'JacobPattern',JStr,'Display',disp,'MaxIter',25,'MaxFunEvals',300*size(xx,2),'TolFun',1e-7,'OutputFcn', @outfun_single_cam);
% options=optimset(options,'JacobPattern',JStr,'Display',disp,'MaxIter',350,'MaxFunEvals',30000*size(xx,2),'TolFun',1e-7);
[chi,resnorm,residuals]=lsqnonlin(@(chi) function_multi_cam_new(xx,chi,X_cell,fc,cc,kc,num_cam),chi0,[],[],options);

% options=optimset('lsqnonlin');
% options=optimset(options,'Display',disp,'MaxIter',0,'MaxFunEvals',30000*size(xx,2),'TolFun',1e-7,'OutputFcn', @outfun_single_cam);
% % % options=optimset(options,'JacobPattern',JStr,'Display',disp,'MaxIter',350,'MaxFunEvals',30000*size(xx,2),'TolFun',1e-7);
% [chi2,resnorm2,residuals2,~,~,~,jacob]=lsqnonlin(@(chi) function_multi_cam_new(xx,chi,X_cell,fc,cc,kc,num_cam),chi0,[],[],options);
% return
% save temp_multi_cam_results chi resnorm residuals
% load temp_multi_cam_results

residuals(residuals==0)=NaN;
residuals_hypot=hypot(residuals(1:2:end,:),residuals(2:2:end,:));
RPE_final=nanmean(residuals_hypot(:));

om_cams_2_ref_opt=chi(:,1:num_cam);
T_cams_2_ref_opt=chi(:,num_cam+1:2*num_cam);

om_grids_2_ref_opt=chi(:,2*num_cam+1:2*num_cam+size(X_cell,2));
T_grids_2_ref_opt=chi(:,2*num_cam+size(X_cell,2)+1:end);

R_ref_cam=rodrigues(om_cams_2_ref_opt(:,ref_cam));
T_ref_cam=T_cams_2_ref_opt(:,ref_cam);
    
om_cams_opt=NaN(3,num_cam);
T_cams_opt=NaN(3,num_cam);
% 
% calib_data.fc=fc;
% calib_data.cc=cc;
% calib_data.kc=kc;
% calib_data.x_cam_all=x_cam;
% calib_data.XX_cam=XX_cam{1};
% calib_data.num_cam=9;
% calib_data.om_cams=om_cams_2_ref_opt;
% calib_data.T_cams=T_cams_2_ref_opt;
% calib_data.om_grids=om_grids_2_ref_opt;
% calib_data.T_grids=T_grids_2_ref_opt;
% calib_data.residuals=residuals;
% calib_data.ref_cam=ref_cam;
% 
% save calib_data_ref calib_data
% return

for c=1:num_cam

    R_cams_2_ref_opt=rodrigues(om_cams_2_ref_opt(:,c));

    R_cams_opt=R_ref_cam'*R_cams_2_ref_opt; % seems correct

    om_cams_opt(:,c)=rodrigues(R_cams_opt);
    T_cams_opt(:,c)=R_ref_cam'*(T_cams_2_ref_opt(:,c)-T_ref_cam); % think is correct

end

clear om_grids T_grids
om_grids_opt=NaN(3,size(om_grids_2_ref_opt,2));
T_grids_opt=NaN(3,size(om_grids_2_ref_opt,2));

for i=1:size(om_grids_2_ref_opt,2)
    
    if ~isnan(om_grids_2_ref_opt(1,i))
        R_grids_2_ref_opt=rodrigues(om_grids_2_ref_opt(:,i));
    
        R_grids_opt=R_ref_cam'*R_grids_2_ref_opt; % seems correct
    
        om_grids_opt(:,i)=rodrigues(R_grids_opt);
        T_grids_opt(:,i)=R_ref_cam'*(T_grids_2_ref_opt(:,i)-T_ref_cam); % think is correct
    end
end

% return

% %% plotting functions for reference
% % plot reference grid
% X_ref=reshape(XX_cam{1}(1,:),nx,[]);
% Y_ref=reshape(XX_cam{1}(2,:),nx,[]);
% Z_ref=reshape(XX_cam{1}(3,:),nx,[]);
% 
% plot3(Z_ref([1:3 nx+1]),X_ref([1:3 nx+1]),Y_ref([1:3 nx+1]),'ok','markersize',3); hold on
% plot3(Z_ref([4:nx nx+2:end]),X_ref([4:nx nx+2:end]),Y_ref([4:nx nx+2:end]),'.k','markersize',1); hold on
% axis image
% 
% working_distance=mean(sqrt(sum(T_cams_2_ref_opt.^2)));
% 
% Cx=[-1 1 1 -1 -1; -0.5 0.5 0.5 -0.5 -0.5]*0.1*working_distance; % makes camera shape for drawing
% Cy=[-1 -1 1 1 -1; -0.5 -0.5 0.5 0.5 -0.5]*0.1*working_distance;
% Cz=[0 0 0 0 0; -1 -1 -1 -1 -1]*0.1*working_distance;
% 
% C=[Cx(:)';Cy(:)';Cz(:)'];
% 
% for c=1:num_cam
%     C2=rodrigues(om_cams_2_ref_opt(:,c))*C+repmat(T_cams_2_ref_opt(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=rand(1,3);
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('Z','X','Y');
%     set(gca,'ydir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% 
% figure
% 
% plot3(Z_ref([1:3 nx+1]),X_ref([1:3 nx+1]),Y_ref([1:3 nx+1]),'ok','markersize',3); hold on
% plot3(Z_ref([4:nx nx+2:end]),X_ref([4:nx nx+2:end]),Y_ref([4:nx nx+2:end]),'.k','markersize',1); hold on
% axis image
% 
% for c=1:num_cam
%     C2=rodrigues(om_cams_2_ref(:,c))*C+repmat(T_cams_2_ref(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=[0 0 0];
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('Z','X','Y');
%     set(gca,'ydir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% return
% 
% for i=1:size(om_grids_2_ref_opt,2)
%     
%     R_grids_2_ref_opt=rodrigues(om_grids_2_ref_opt(:,i));
%     
%     XX2=R_grids_2_ref_opt*XX_cam{i}+repmat(T_grids_2_ref_opt(:,i),1,size(XX_cam{i},2));
% 
%     X2=reshape(XX2(1,:),nx,[]);
%     Y2=reshape(XX2(2,:),nx,[]);
%     Z2=reshape(XX2(3,:),nx,[]);
% 
%     col=rand(1,3);
% 
%     p1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o','color',col); hold on
%     p2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.','color',col); hold on
% %             
% end

% return

% figure
% for c=1:num_cam
%     C2=rodrigues(om_cams_opt(:,c))*C+repmat(T_cams_opt(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=rand(1,3);
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('Z','X','Y');
%     set(gca,'ydir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% 
% for i=1:size(om_grids_2_ref_opt,2)
%     
%     R_grids_opt=rodrigues(om_grids_opt(:,i));
%     
%     XX2=R_grids_opt*XX_cam{i}+repmat(T_grids_opt(:,i),1,size(XX_cam{i},2));
% 
%     X2=reshape(XX2(1,:),nx,[]);
%     Y2=reshape(XX2(2,:),nx,[]);
%     Z2=reshape(XX2(3,:),nx,[]);
% 
%     col=rand(1,3);
% 
%     p1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o','color',col); hold on
%     p2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.','color',col); hold on
% %             
% end