function [fc_opt, cc_opt, T_opt, om_opt, residual]=single_cam_calib_cc(fc,cc,kc,T,om,x_cam,XX_cam,disp)
% save temp_single_cam2
% moo=1
% return
% clear 
% close all force
% load temp_single_cam2

% clear disp
% disp='off';

% use=[1:2];
% 
% x_cam=x_cam(use);
% XX_cam=XX_cam(use);
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

xx=cell2mat(x_cell);
grid_t=numel(xx);

%BUNDLE ADJUSTMENT ROUTINE
%Specifies structure of the Jacobian matrix in the sparse matrix JStr

clear j*
chi0=[[fc(1);cc(1);cc(2)] om T];
jcol_cam=repmat([1 2 1 3],size(xx,2),1);

jcol_cam=jcol_cam(:);
% jrow_cam=[1:2:2*size(xx,2) 2:2:2*size(xx,2)]';
% jrow_cam=[1:2:2*size(xx,2) 2:2:2*size(xx,2) 1:2:2*size(xx,2)]';
jrow_cam=[1:2:2*size(xx,2) 1:2:2*size(xx,2) 2:2:2*size(xx,2) 2:2:2*size(xx,2)]';

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

% return

% %Nonlinear least squares optimization routine for bundle adjustment
options=optimset('lsqnonlin');
options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam);
[chi,resnorm,residual,exitflag,output]=lsqnonlin(@(chi) function_single_cam_cc(xx,chi,X_cell,kc),chi0,[],[],options);

% options=optimset('lsqnonlin');
% options=optimset(options,'Display','off','MaxIter',1,'MaxFunEvals',200*length(xx));%,'OutputFcn', @outfun_single_cam);
% [chi2,resnorm2,residual,exitflag,output,lambda,jacob]=lsqnonlin(@(chi) function_single_cam_cc(xx,chi,X_cell,kc),chi0,[],[],options);

fc_opt=chi([1 1],1);
cc_opt=chi([2 3],1);

om_opt=chi(:,2:size(X_cell,2)+1);
T_opt=chi(:,size(X_cell,2)+2:end);

%     disp(resnorm)
%     return
