function [cam_distance, cam_angle]=plot_calibration_setup_function(calib_data)
% save temp_plot_calib
% clear
% close all
% load temp_plot_calib

num_cam=size(calib_data.om_cams,2);

om_cams_opt=calib_data.om_cams;
T_cams_opt=calib_data.T_cams;

om_grids_opt=calib_data.om_grids;
T_grids_opt=calib_data.T_grids;
XX_cam=calib_data.XX_cam;

R2=rodrigues(om_grids_opt(:,1));
T2=T_grids_opt(:,1);

R21=R2';
T21=-R21*T2; 

orientate_2_grid=0;

if isfield(calib_data,'ref_grid')
    ref_grid=calib_data.ref_grid;
else
    ref_grid=1;
end

nx=calib_data.nx;

% nx=sqrt(size(XX_cam,2));

T_grids_dist=nanmean(T_grids_opt(3,:));

working_distance=nanmean(T_grids_opt(3,:));

Cx=[-1 1 1 -1 -1; -0.5 0.5 0.5 -0.5 -0.5]*0.1*working_distance; % makes camera shape for drawing
Cy=[1 1 -1 -1 1; 0.5 0.5 -0.5 -0.5 0.5]*0.1*working_distance;
Cz=[0 0 0 0 0; -1 -1 -1 -1 -1]*0.1*working_distance;

if isfield(calib_data,'ref_cam')
   ref_cam=calib_data.ref_cam;
else
    ref_cam=1;
end
    
figure
for i=1:num_cam
    
    T1=T_cams_opt(:,i);
    R=rodrigues(om_cams_opt(:,i));
    C=[Cx(:)';Cy(:)';Cz(:)'];
    C2=R*C+repmat(T1,1,size(C,2));
    
    C0=R*[0;0;0]+T1;
    
    
    if orientate_2_grid==1
        C2=R21*C2+repmat(T21,1,size(C2,2));
        T_2_cent(:,i)=R21*C0+T21;
    end
    
    Cx2=reshape(C2(1,:),2,[]);
    Cy2=reshape(C2(2,:),2,[]);
    Cz2=reshape(C2(3,:),2,[]);
    
    Cx2_mean=mean(Cx2(:,1:4),2);
    Cy2_mean=mean(Cy2(:,1:4),2);
    Cz2_mean=mean(Cz2(:,1:4),2);
   
    C_vec1(i,:)=[Cx2_mean(1) Cy2_mean(1) Cz2_mean(1)];
    C_vec2(i,:)=[Cx2_mean(2) Cy2_mean(2) Cz2_mean(2)];
   
    if orientate_2_grid==1
        text_pos(i,:)=[mean(Cx2(:)),mean(Cy2(:)),mean(Cz2(:))+0.15*working_distance];
        text_pos2(i,:)=[mean(Cx2(:)),mean(Cy2(:)),mean(Cz2(:))+0.21*working_distance];

        h1=plot3(Cx2,Cy2,Cz2,Cx2',Cy2',Cz2'); hold on; axis image
        p=patch(Cx2([5 6 8 7]),Cy2([5 6 8 7]),Cz2([5 6 8 7]),'b','edgecolor','none'); axis image
    else
        text_pos(i,:)=[mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance];
        text_pos2(i,:)=[mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.21*working_distance];

        h1=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
        p=patch(Cz2([5 6 8 7]),Cx2([5 6 8 7]),Cy2([5 6 8 7]),'b','edgecolor','none'); axis image
    end
    
    if i==ref_cam
        set(h1,'color','r','linewidth',2)
        set(p,'facecolor','r')
    else
       set(h1,'color','b','linewidth',2);
       set(p,'facecolor','b')
    end
   
end

for i=1:size(om_grids_opt,2)
    XX2=rodrigues(om_grids_opt(:,i))*XX_cam+repmat(T_grids_opt(:,i),1,size(XX_cam,2));
    
    if orientate_2_grid==1
        XX2=R21*XX2+repmat(T21,1,size(XX2,2));
    end
    
    X2=reshape(XX2(1,:),nx,[]);
    Y2=reshape(XX2(2,:),nx,[]);
    Z2=reshape(XX2(3,:),nx,[]);
 
    col=rand(1,3);
    
    if orientate_2_grid==1
        h1=plot3(X2([1:3 nx+1]),Y2([1:3 nx+1]),Z2([1:3 nx+1]),'o','color',col,'markersize',4); hold on
        h2=plot3(X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),Z2([4:nx nx+2:end]),'.','color',col,'markersize',5); hold on
    else
        h1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o','color',col,'markersize',4); hold on
        h2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.','color',col,'markersize',5); hold on
    end
end

if orientate_2_grid==1
    title('Reference camera axis system')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
else
    title('Reference camera axis system')
    xlabel('Z (mm)')
    ylabel('X (mm)')
    zlabel('Y (mm)')
    
    set(gca,'ydir','reverse','zdir','reverse')

end


axis image
% close all

[P_intersect,~] = lineIntersect3D(C_vec1,C_vec2);

% plot3(P_intersect(3),P_intersect(1),P_intersect(2),'.r','markersize',20)
% plot3(P_intersect([3 3]),[0 0],zlim,'r');
% plot3(P_intersect([3 3]),ylim,[0 0],'r');
% plot3(xlim,[0 0],[0 0],'r');
  
if orientate_2_grid==1
    cam_distance=sqrt(sum(T_2_cent.^2));
    cam_angle=atan2d(T_2_cent(3,:),hypot(T_2_cent(1,:),T_2_cent(2,:)));
else
    T_2_cent=T_cams_opt-repmat([0;0;P_intersect(3)],1,num_cam);

    cam_distance=sqrt(sum(T_2_cent.^2));
    cam_angle=atan2d(-T_2_cent(2,:),hypot(T_2_cent(1,:),T_2_cent(3,:)));
end

cam_dist2=round(cam_distance);
cam_angle2=0.1*round(10.*cam_angle);

for i=1:num_cam
   
    t1=text(text_pos(i,1),text_pos(i,2),text_pos(i,3),['Camera ' num2str(i)]);
    t2=text(text_pos2(i,1),text_pos2(i,2),text_pos2(i,3),['(' num2str(cam_angle2(i)) '^o, ' num2str(cam_dist2(i)) 'mm)']);
    
    if orientate_2_grid==1
        p1=plot3([C_vec1(i,1) P_intersect(3)],...
                 [C_vec1(i,2) 0],...
                 [C_vec1(i,3) 0],'--r');
    else
        p1=plot3([C_vec1(i,3) P_intersect(3)],...
                 [C_vec1(i,1) 0],...
                 [C_vec1(i,2) 0],'--r');
    end
    
      
    if i==ref_cam
        set(t1,'color','r')
        set(t2,'color','r')
        set(p1,'color','r')
    else
        set(t1,'color','b')
        set(t2,'color','b')
        set(p1,'color','b')
    end
end