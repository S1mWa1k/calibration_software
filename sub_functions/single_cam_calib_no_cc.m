function [fc_opt, T_opt, om_opt, x_cam, removed, residual]=single_cam_calib_no_cc(fc,cc,kc,T,om,x_cam,XX_cam)
% save temp_single_cam
% moo=1
% return
% clear 
% close all force
% load temp_single_cam

% clear disp
% disp='off';

% use=[1:3];
% x_cam=x_cam(use);
% XX_cam=XX_cam(use);
% 
% x_cam{1}(:,3:5)=NaN;
% x_cam{2}(:,8)=NaN;
% 
% T=T(:,use);
% om=om(:,use);

k=1;
for i=1:length(x_cam)
    
    xx=x_cam{i};
%     xx=xx(:,1:10);
    
    if ~isempty(xx)
        [Ay, Ax]=find(isnan(xx));
        xx(:,Ax)=[];

        
        X=XX_cam{i};
%         X=X(:,1:10);
        X(:,Ax)=[];
        
        x_cell{k}=xx;
        X_cell{k}=X;
        
        grid_l(i)=size(X,2);
        
        k=k+1;
    end
    
end

%BUNDLE ADJUSTMENT ROUTINE
%Specifies structure of the Jacobian matrix in the sparse matrix JStr

changed=0;
removed=zeros(size(x_cam));

for runs=1:2

    if runs==1 || (runs==2 && changed==1)
    
        xx=cell2mat(x_cell);
        grid_t=numel(xx);
    
        clear j*
     
        chi0=[[fc(1);0;0] om T];
        jcol_cam=repmat([1 1],size(xx,2),1);
    
        jcol_cam=jcol_cam(:);
        jrow_cam=[1:2:2*size(xx,2) 2:2:2*size(xx,2)]';
    
        for i=1:length(grid_l)
            jcol_om_grid{i,1}=repmat(3+(1:3)+3*(i-1),2*grid_l(i),1);
            jrow_om_grid{i,1}=repmat((1:2*grid_l(i))+2*sum(grid_l(1:i-1)),3,1)';
    
            jcol_T_grid{i,1}=jcol_om_grid{i,1}+3*length(grid_l);
            jcol_T_grid{i,1}(2:2:size(jcol_T_grid{i,1},1),1)=NaN;
            jcol_T_grid{i,1}(1:2:size(jcol_T_grid{i,1},1),2)=NaN;
    
            jrow_T_grid{i,1}=jrow_om_grid{i,1};
            jrow_T_grid{i,1}(2:2:size(jcol_T_grid{i,1},1),1)=NaN;
            jrow_T_grid{i,1}(1:2:size(jcol_T_grid{i,1},1),2)=NaN;
    
        end
    
        jcol_om_grid=cell2mat(jcol_om_grid);
        jrow_om_grid=cell2mat(jrow_om_grid);
    
        jcol_T_grid=cell2mat(jcol_T_grid);
        jrow_T_grid=cell2mat(jrow_T_grid);
    
        jcol_T_grid(isnan(jcol_T_grid))=[];
        jrow_T_grid(isnan(jrow_T_grid))=[];
    
        jcol=[jcol_cam(:);jcol_om_grid(:);jcol_T_grid(:)];
        jrow=[jrow_cam(:);jrow_om_grid(:);jrow_T_grid(:)];
    
        JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,3+6*length(grid_l));
    
        % %Nonlinear least squares optimization routine for bundle adjustment
        options=optimset('lsqnonlin');
        options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam);
        [chi,resnorm,residual,exitflag,output]=lsqnonlin(@(chi) function_single_cam_no_cc(xx,chi,X_cell,cc,kc),chi0,[],[],options);
    
    %     options=optimset('lsqnonlin');
    %     options=optimset(options,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx));%,'OutputFcn', @outfun_single_cam);
    %     [chi,resnorm,residual,exitflag,output,lambda,jacob]=lsqnonlin(@(chi) function_single_cam_no_cc(xx,chi,X_cell,cc,kc,runs),chi0,[],[],options);
        
    %     disp(resnorm)
        
        if runs==1
            for i=1:length(x_cell)
        
                if ~isempty(x_cell{i})
            
                    xx_temp=x_cell{i};
            
                    F=function_single_cam_no_cc(xx_temp,chi(:,[1 i+1 length(x_cell)+i+1]),X_cell(i),cc,kc);
        
                    RR=hypot(F(1,:),F(2,:));
                    
                    MSE=mean(RR);
        
                    if MSE>5
        
                        changed=1;
        
                        x_cell{i}=[];
                        X_cell{i}=double.empty(3,0);

                        om(:,i)=NaN(3,1);
                        T(:,i)=NaN(3,1);

                        grid_l(i)=0;

                        x_cam{i}=NaN(size(x_cam{i}));

                        removed(i)=1;
        
                    end
        
                    
                end
        
            end
        end
    end

end

fc_opt=chi([1 1],1);

om_opt=chi(:,2:size(X_cell,2)+1);
T_opt=chi(:,size(X_cell,2)+2:end);

RR=hypot(residual(1,:),residual(2,:));
MSE=mean(RR);