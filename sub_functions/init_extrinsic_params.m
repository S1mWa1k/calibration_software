function [om, T,x_cam1, removed, residuals_MSE, residuals_MaxSE]=init_extrinsic_params(x_cam1,XX_cam,fc,cc,kc,alpha_c)
% save temp_init_extrinsic

% clear all
% close all
% load temp_init_extrinsic
% clc

%%% Computes the extrinsic parameters for all the active calibration images 

residuals_MSE=NaN(1,length(x_cam1));
residuals_MaxSE=NaN(1,length(x_cam1));
removed=zeros(1,length(x_cam1));

for kk = 1:length(x_cam1)

    if any(~isnan(x_cam1{kk}(1,:)))
        changed=0;
    
        for run=1:2
    
            if run==1 || (run==2 && changed==1)
                
                x_kk=x_cam1{kk};
                X_kk=XX_cam{kk};
                
                [r c]=find(isnan(x_kk));
                x_kk(:,c)=[];
                X_kk(:,c)=[];
                
                [r c]=find(isnan(X_kk));
                x_kk(:,c)=[];
                X_kk(:,c)=[];    
            %     save temp
                if isempty(x_kk)
                    omckk=NaN(3,1);
                    Tckk=NaN(3,1);
                else
                    x_kk=x_kk(1:2,:);
                    X_kk=[X_kk(1:2,:);ones(1,size(X_kk,2))];
                    
                    N_points_views(kk) = size(x_kk,2);
                    [omckk,Tckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c);
                    [omckk,Tckk,Rckk,JJ_kk] = compute_extrinsic_refine(omckk,Tckk,x_kk,X_kk,fc,cc,kc,alpha_c,20,1e6);
                end
            %     save temp
                
                if~isnan(omckk)
                    Rckk=rodrigues(omckk);
                    om_temp=rodrigues(Rckk);
            %         Rckk(1:2,:)=-Rckk(1:2,:);
                    om(:,kk)=rodrigues(Rckk);
                else
                    om_temp=omckk;
                    om(:,kk)=omckk;
                end
                
                Tckk_temp=Tckk;
            %     Tckk(1:2)=-Tckk(1:2);
                T(:,kk)=Tckk;
                
                xp = project_points2(X_kk,om_temp,Tckk_temp,fc,cc,kc,alpha_c);
            
                residuals_temp=[xp(1,:)-x_kk(1,:);xp(2,:)-x_kk(2,:)];
            
                RR=hypot(residuals_temp(1:2:end,:),residuals_temp(2:2:end,:));
            
                % residuals_temp{kk}=[xp(1,:)-x_kk(1,:);xp(2,:)-x_kk(2,:)];
        
                if ~isempty(RR) && run==1
                          
                    if any(RR>20)
    
                        x_cam1{kk}(:,RR>20)=NaN;
                        changed=1;
                   
                    end
                                
                end
            end
        end
    
        if sum(RR>10)>10 || sum(~isnan(x_cam1{kk}(1,:)))<30
            
            x_cam1{kk}=NaN(size(x_cam1{kk}));
        
            om(:,kk)=NaN(3,1);
            T(:,kk)=NaN(3,1);
        
            removed(kk)=1;
        
        else
        
            residuals_MSE(1,kk)=mean(RR,'omitnan');
            residuals_MaxSE(1,kk)=max(RR);
        
            % if changed==1
            %     disp([residuals_MSE_init residuals_MSE(kk)])
            % end
        
        end
    else

        residuals_MSE(1,kk)=NaN;
        residuals_MaxSE(1,kk)=NaN;
        om(:,kk)=NaN(3,1);
        T(:,kk)=NaN(3,1);

    end

end

% residuals=residuals_temp;
% residuals=cell2mat(residuals_temp);