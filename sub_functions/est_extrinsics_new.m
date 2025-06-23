function [T_cams_2_ref, om_cams_2_ref, om_grids_2_ref, T_grids_2_ref, ref_grid, xx_cam]=est_extrinsics_new(fc,cc,kc,om,T,xx_cam,XX_cam,draw)
% save temp_est_extrinsics
% save temp_est_extrinsics_0dist

% clear
% close all
% load temp_est_extrinsics
% load temp_est_extrinsics_0dist

% % % 
% draw='off';
% nx=14;

n_cams=length(om);
n_grids=length(XX_cam);

grid_count=zeros(n_cams,n_grids);
% return
working_distance=900;

Cx=[-1 1 1 -1 -1; -0.5 0.5 0.5 -0.5 -0.5]*0.1*working_distance; % makes camera shape for drawing
Cy=[-1 -1 1 1 -1; -0.5 -0.5 0.5 0.5 -0.5]*0.1*working_distance;
Cz=[0 0 0 0 0; -1 -1 -1 -1 -1]*0.1*working_distance;

C=[Cx(:)';Cy(:)';Cz(:)'];

% First find if there are any grids that are seen in all camera views.
% Also reshape om and T so that they are grouped into grid number

residuals_ssq=NaN(n_cams,n_grids);
n_res=NaN(n_cams,n_grids);

for c=1:n_cams

    grid_count(c,:)=~isnan(om{c}(1,:));
    
    cx=fc(1,c);
    cy=fc(2,c);

    K=kc(1:3,c);
    P=kc(4:5,c);
        
    for j=1:n_grids
       
        om_array{c,j}=om{c}(:,j);
        T_array{c,j}=T{c}(:,j);
        
        R2=rodrigues(om{c}(:,j));
        T2=T{c}(:,j);
       
        XWP=-R2'*T2-XX_cam{j};
        
        XWP2=R2(3,:)*XWP;

        if ~isnan(T{c}(1,j))

            dx_e=cx*(R2(1,:)*XWP)./XWP2;
            dy_e=cy*(R2(2,:)*XWP)./XWP2;

            dx_ne=dx_e/cx;
            dy_ne=dy_e/cy;

            r_ne=sqrt(dx_ne.^2+dy_ne.^2);

            n=size(XWP,2);
            
            im_x=cc(1,c)+dx_e./(ones(1,n)-(K(1).*r_ne.^2+K(2).*r_ne.^4+K(3).*r_ne.^6)+2*P(1).*dx_ne.*dy_ne+P(2)*(r_ne.^2+2*dx_ne.^2));
            im_y=cc(2,c)+dy_e./(ones(1,n)-(K(1).*r_ne.^2+K(2).*r_ne.^4+K(3).*r_ne.^6)+2*P(2).*dx_ne.*dy_ne+P(1)*(r_ne.^2+2*dy_ne.^2));

            % im_x=cx*(R2(1,:)*XWP)./XWP2+cc(1,c);
            % im_y=cy*(R2(2,:)*XWP)./XWP2+cc(2,c);
    
            residuals=[im_x-xx_cam{c,j}(1,:); im_y-xx_cam{c,j}(2,:)];
            
            RR=hypot(residuals(1,:),residuals(2,:));
            
            xx_cam{c,j}(:,RR>3)=NaN;

            if sum(~isnan(xx_cam{c,j}(1,:)))<30

                xx_cam{c,j}=NaN(size(xx_cam{c,j}));

                om{c}(:,j)=NaN(3,1);
                T{c}(:,j)=NaN(3,1);

                om_array{c,j}=NaN(3,1);
                T_array{c,j}=NaN(3,1);

                RR=NaN(size(RR));
                
            end

            residuals_ssq(c,j)=nansum(RR.^2);

            % n_res(c,j)=sum(~isnan(residuals));
            
            % plot(im_x,im_y,'ow'); hold on
            % plot(xx_cam{c,j}(1,:),xx_cam{c,j}(2,:),'.g','markersize',12);
            % 
            % plot([im_x;xx_cam{c,j}(1,:)],[im_y;xx_cam{c,j}(2,:)],'linewidth',3,'color',[0.8 0 0])
            % set(gca,'color','k')
            % axis image
            % return
        end
        
    end
   
end

res_var=nansum(residuals_ssq)./(nansum(n_res)-1);
% return
grid_score=sum(grid_count);

% find a grid that contrains all cameras views (assuming there is
% one). If more than one, find one with lowest residual variance
pot_ref_grid=find(grid_score==max(grid_score));

ref_grid=pot_ref_grid(find(res_var(pot_ref_grid)==min(res_var(pot_ref_grid)),1));
% ref_grid=1;
% can now work out the orientation and tranlations of the other grids
% relative to the reference grid
om_ref_grid=reshape(cell2mat(om_array(:,ref_grid)),3,[]);
T_ref_grid=reshape(cell2mat(T_array(:,ref_grid)),3,[]);

% % plot grids for one camera, just for viewing pleasure
% col=[0.7 0.2 0;
%      0.0 0.2 0.7];
% 
% for c=1
% 
%     if ~isnan(om{c}(1,1)) && ~isnan(om{c}(1,ref_grid))
% 
%         figure
% 
%         for i=[1  ref_grid]
% 
%             XX2=rodrigues(om{c}(:,i))*XX_cam{i}+repmat(T{c}(:,i),1,size(XX_cam{i},2));
% 
%             X2=reshape(XX2(1,:),nx,[]);
%             Y2=reshape(XX2(2,:),nx,[]);
%             Z2=reshape(XX2(3,:),nx,[]);
% 
%             p1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o'); hold on
%             p2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.'); hold on
% 
%             if i~=ref_grid
%                 set(p1,'color',[0.8 0.1 0])
%                 set(p2,'color',[0.8 0.1 0])
%             elseif i==ref_grid
%                 set(p1,'color','k')
%                 set(p2,'color','k')
%             end
%         end
% 
%         Cx=reshape(C(1,:),2,[]);
%         Cy=reshape(C(2,:),2,[]);
%         Cz=reshape(C(3,:),2,[]);
% 
%         h=plot3(Cz,Cx,Cy,Cz',Cx',Cy'); hold on; axis image
%         set(h,'color','b','linewidth',2);
%         h2=text(mean(Cz(:)),mean(Cx(:)),mean(Cy(:))+0.15*working_distance,num2str(i));
%         set(h2,'color','b')
%         patch(Cz([1 2 4 3]),Cx([1 2 4 3]),Cy([1 2 4 3]),'b'); axis image
% 
%         plot3(mean(Cz(:,1:4),2),mean(Cx(:,1:4),2),mean(Cy(:,1:4),2),'b'); hold on; axis image
% 
%         plot3(0,0,0,'.k');
% %         view(-163,-90)
%         labels('Z','X','Y');
%         set(gca,'ydir','reverse','zdir','reverse')
% 
%         axis image
%     end
% end
% return
% figure

om_grids_2_ref=NaN(3,n_grids);
T_grids_2_ref=NaN(3,n_grids);

for i=1:n_grids
   
    if i==ref_grid
        
        om_grids_2_ref(:,i)=0;
        T_grids_2_ref(:,i)=0;

        for i=1:n_cams
            om_grid_2_ref_all{i,c}=NaN(3,1);
            T_grid_2_ref_all{i,c}=NaN(3,1);
        end
    else
        
        om_temp=reshape(cell2mat(om_array(:,i)),3,[]);
        T_temp=reshape(cell2mat(T_array(:,i)),3,[]);
                
        om_grid_2_ref_temp=NaN(3,n_cams);
        T_grid_2_ref_temp=NaN(3,n_cams);
%         close
        for c=1:n_cams
            
            if ~isnan(T_temp(1,c)) && ~isnan(T_ref_grid(1,c))

                R_ref_grid=rodrigues(om_ref_grid(:,c));
                R_temp=rodrigues(om_temp(:,c));

                R_grid_2_ref=R_ref_grid'*R_temp; % think is correct

                om_grid_2_ref_temp(:,c)=rodrigues(R_grid_2_ref);
                T_grid_2_ref_temp(:,c)=R_ref_grid'*(T_temp(:,c)-T_ref_grid(:,c)); % think is correct
                
                om_grid_2_ref_all{i,c}=om_grid_2_ref_temp;
                T_grid_2_ref_all{i,c}=T_grid_2_ref_temp;

                XX2=R_grid_2_ref*XX_cam{i}+repmat(T_grid_2_ref_temp(:,c),1,size(XX_cam{i},2));

%                 X2=reshape(XX2(1,:),nx,[]);
%                 Y2=reshape(XX2(2,:),nx,[]);
%                 Z2=reshape(XX2(3,:),nx,[]);
%                 
%                 col=rand(1,3);

%                 p1=plot3(X2([1:3 nx+1]),Y2([1:3 nx+1]),Z2([1:3 nx+1]),'o','color',col); hold on
%                 p2=plot3(X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),Z2([4:nx nx+2:end]),'.','color',col); hold on
%             
            else
                
                om_grid_2_ref_all{i,c}=NaN(3,1);
                T_grid_2_ref_all{i,c}=NaN(3,1);
            end
                        
%             X_ref=reshape(XX_cam{i}(1,:),nx,[]);
%             Y_ref=reshape(XX_cam{i}(2,:),nx,[]);
%             Z_ref=reshape(XX_cam{i}(3,:),nx,[]);
%             
%             plot3(X_ref([1:3 nx+1]),Y_ref([1:3 nx+1]),Z_ref([1:3 nx+1]),'ok','markersize',5); hold on
%             plot3(X_ref([4:nx nx+2:end]),Y_ref([4:nx nx+2:end]),Z_ref([4:nx nx+2:end]),'.k','markersize',10); hold on
%         
%             axis image
%             view(90,90)

        end
        
        om_grids_2_ref(:,i)=nanmedian(om_grid_2_ref_temp,2);
        T_grids_2_ref(:,i)=nanmedian(T_grid_2_ref_temp,2);
        
%         set(gca,'ydir','reverse','zdir','reverse')
        
    end
    
end

om_cams_2_ref=NaN(3,n_cams);
T_cams_2_ref=NaN(3,n_cams);

if strcmp(draw,'on')
    figure
end

for c=1:n_cams
    
    om_temp=reshape(cell2mat(om_array(c,:)),3,[]);
    T_temp=reshape(cell2mat(T_array(c,:)),3,[]);
    
    om_cam_2_ref_temp=NaN(3,n_grids);
    T_cam_2_ref_temp=NaN(3,n_grids);
    
%     subplot(3,3,c)
    for i=1:n_grids
    
        if ~isnan(om_temp(1,i)) && ~isnan(om_grids_2_ref(1,i))
            R_grid_2_ref=rodrigues(om_grids_2_ref(:,i));
            R_temp=rodrigues(om_temp(:,i));

            % reverse of grid to camera
            R_cam_2_ref_temp=R_grid_2_ref*R_temp';
            om_cam_2_ref_temp(:,i)=rodrigues(R_cam_2_ref_temp);
            T_cam_2_ref_temp(:,i)=-R_cam_2_ref_temp*T_temp(:,i)+T_grids_2_ref(:,i);
            
%             % plot cameras to check they're all in similar positions
%             C2=R_cam_2_ref_temp*C+repmat(T_cam_2_ref_temp(:,i),1,size(C,2));
%             
%             Cx2=reshape(C2(1,:),2,[]);
%             Cy2=reshape(C2(2,:),2,[]);
%             Cz2=reshape(C2(3,:),2,[]);
%             
%             col=rand(1,3);
%             
%             h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%             set(h,'color',col,'linewidth',2);
%             h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(i));
%             set(h2,'color',col)
% %             patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
%             patch(Cz2([5 6 8 7]),Cx2([5 6 8 7]),Cy2([5 6 8 7]),col); axis image
%             
%             h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%             set(h,'color',col,'linewidth',2);
            
        end
    end
    
%     return
    om_cams_2_ref(:,c)=nanmedian(om_cam_2_ref_temp,2);
    T_cams_2_ref(:,c)=nanmedian(T_cam_2_ref_temp,2);
    
    if strcmp(draw,'on')
        % plot reference grid
        X_ref=reshape(XX_cam{i}(1,:),nx,[]);
        Y_ref=reshape(XX_cam{i}(2,:),nx,[]);
        Z_ref=reshape(XX_cam{i}(3,:),nx,[]);

        plot3(Z_ref([1:3 nx+1]),X_ref([1:3 nx+1]),Y_ref([1:3 nx+1]),'ok','markersize',5); hold on
        plot3(Z_ref([4:nx nx+2:end]),X_ref([4:nx nx+2:end]),Y_ref([4:nx nx+2:end]),'.k','markersize',10); hold on
        axis image

        C2=rodrigues(om_cams_2_ref(:,c))*C+repmat(T_cams_2_ref(:,c),1,size(C,2));

        Cx2=reshape(C2(1,:),2,[]);
        Cy2=reshape(C2(2,:),2,[]);
        Cz2=reshape(C2(3,:),2,[]);

        col=rand(1,3);

        h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
        set(h,'color','k','linewidth',2);
        h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
        set(h2,'color','k')
%         patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
        patch(Cz2([5 6 8 7]),Cx2([5 6 8 7]),Cy2([5 6 8 7]),col); axis image

        h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
        set(h,'color','k','linewidth',2);

        labels('Z','X','Y');
        set(gca,'ydir','reverse','zdir','reverse')

        axis image
        labels
    end
%     return
end 
