function [fc_opt, cc_opt, kc_opt, T_opt, om_opt, x_cam, removed, residual]=single_cam_calib_cc_dist(fc,cc,kc,T,om,x_cam,XX_cam)

% save temp_single_cam

% clear 
% close all force
% load temp_single_cam

% clear disp
% disp='off';

% use=1;
% 
% x_cam=x_cam(use);
% XX_cam=XX_cam(use);
% 
% T=T(:,use);
% om=om(:,use);

% then with distortion
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

changed=0;
removed=zeros(size(x_cam));

for runs=1:2

    if runs==1 || (runs==2 && changed==1)

        xx=cell2mat(x_cell);
        grid_t=numel(xx);
        
        %BUNDLE ADJUSTMENT ROUTINE
        %Specifies structure of the Jacobian matrix in the sparse matrix JStr
        
        clear j*
        chi0=[[fc;0] [cc;0] zeros(3,1) om T];
        % chi0=[[fc(1); 0;0] [cc;0] zeros(3,1) om T];

        % jcol_cam=repmat([1 2 4 5 7 8 1 2 4 5 7 8],size(xx,2),1);
        jcol_cam=repmat([1 4 5 7 8 1 4 5 7 8],size(xx,2),1);
        
        jrow_cam=[repmat(1:2:2*size(xx,2),1,0.5*size(jcol_cam,2)) repmat(2:2:2*size(xx,2),1,0.5*size(jcol_cam,2))]';
        jcol_cam=jcol_cam(:);
        
        % return
        % 2:2:2*size(xx,2)
        for i=1:length(grid_l)
            jcol_om_grid{i,1}=repmat(9+(1:3)+3*(i-1),2*grid_l(i),1);
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
        
        JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,9+6*length(grid_l));
        
        % %Nonlinear least squares optimization routine for bundle adjustment
        options=optimset('lsqnonlin');
        options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam);
        [chi,resnorm,residual,exitflag,output]=lsqnonlin(@(chi) function_single_cam_cc_dist(xx,chi,X_cell),chi0,[],[],options);
        
        % options=optimset('lsqnonlin');
        % options=optimset(options,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam);%,'OutputFcn', @outfun_single_cam);
        % [chi2,resnorm2,residual2,exitflag,output,lambda,jacob]=lsqnonlin(@(chi) function_single_cam_cc_dist(xx,chi,X_cell),chi0,[],[],options);
        % % 
        % spy(JStr,'r',25); axis equal
        % hold on
        % figure;
        % spy(jacob,'b',11); axis equal
        % return

        if runs==1
            for i=1:length(x_cell)
        
                if ~isempty(x_cell{i})
            
                    xx_temp=x_cell{i};
            
                    F=function_single_cam_cc_dist(xx_temp,chi(:,[1:3 i+3 length(x_cell)+i+3]),X_cell(i));
        
                    RR=hypot(F(1,:),F(2,:));
                    
                    MSE=mean(RR);
        
                    if MSE>3
        
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
% fc_opt=chi([1 2],1);
cc_opt=chi([1 2],2);

kc_opt=[chi(1:2,3);0;0;0];

om_opt=chi(:,4:size(X_cell,2)+3);
T_opt=chi(:,size(X_cell,2)+4:end);
