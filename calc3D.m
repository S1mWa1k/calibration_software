function [X,Y,Z,resnorm,residual_all,coords,X0,Y0,Z0,x_rep,y_rep]=calc3D(xx,yy,calib_data)

% INPUTS
% xx, yy are the 'x' and 'y' coordinates, these must be in a nxm matrix
% where n is the number of cameras, and m is the number of datapoints
%
% calib_data is a calibration data file
%
% num_cam is the number of cameras used
%
% OUTPUTS
% X, Y, Z, are the 3D coordinates, after rotating the data according the
% orientation and position of the first grid, which is normally aligned with gravity
%
% resnorm is the sum of squares of the residuals
%
% residual_all are the residuals
%
% coords are just the same as the input but reordered into a single
% parameter
%
% X0, Y0, Z0, are the 3D coordinates, but in the coordinate system of the
% first camera
%
% x_rep and y_rep are the reprojected coordinates

num_cam=size(xx,1);

coords=NaN(2*num_cam,size(xx,2));

for i=1:num_cam
    
    coords(2*i-1,:)=xx(i,:);
    coords(2*i,:)=yy(i,:);
    
end

x=coords;

fc=calib_data.fc;
cc=calib_data.cc;

T_cams=calib_data.T_cams;
om_cams=calib_data.om_cams;

for i=1:num_cam
    K{i}=[fc(1,i) 0 cc(1,i);
        0 fc(2,i) cc(2,i);
        0 0 1];

    R{i}=rodrigues(om_cams(:,i))';
    T_cams2(:,i)=R{i}*T_cams(:,i);

    P{i}=K{i}*[R{i} T_cams2(:,i)];
end

XX=zeros(3,size(x,2));

for i=1:size(x,2)
    
    for j=1:num_cam
        
        A(2*j-1,:)=x(2*j-1,i)*P{j}(3,:)-P{j}(1,:);
        A(2*j,:)=x(2*j,i)*P{j}(3,:)-P{j}(2,:);
        
    end
       
   A(isnan(A(:,1)),:)=[];
   
   if size(A,1)>2
       [~, ~, v]=svd(A);
       XX(:,i)=v(1:3,end)./-v(4,end);
   else
       XX(:,i)=NaN;
   end
end

XX(XX==0)=NaN;

for c=1:num_cam
    c_x=fc(1,c);
    c_y=fc(2,c);                        
    
    T=repmat(T_cams(:,c),1,size(XX,2));
    R=rodrigues(om_cams(:,c))';  
    xr(2*c-1,:)=c_x*(R(1,:)*(T-XX))./(R(3,:)*(T-XX))+cc(1,c);
    xr(2*c,:)=c_y*(R(2,:)*(T-XX))./(R(3,:)*(T-XX))+cc(2,c);
    
%     subplot(2,2,c); hold on
%     plot(xx(c,:),yy(c,:),'.'); axis image
%     plot(xr(2*c-1,:),xr(2*c,:),'.g'); axis image
end

x_rep=xr(1:2:end,:);
y_rep=xr(2:2:end,:);

residual_all=xr-coords;
resnorm=nansum(nansum(residual_all.^2));

if isfield(calib_data,'om_grids')

    R2=rodrigues(calib_data.om_grids(:,1));
    T2=calib_data.T_grids(:,1);
    rotate=0;
else
    R2=rodrigues(calib_data.om_bar(:,1));
    T2=calib_data.T_bar(:,1);
    rotate=1;
end

R1=rodrigues([0;0;0]);
T1=[0;0;0];

R21=R1*R2';
T21=-R1*R2'*T2+T1;

XX2=R21*XX+repmat(T21,1,size(XX,2));

X=XX2(1,:);
Y=XX2(2,:);
Z=XX2(3,:);

if rotate==1
   [X,Y,Z]= Ry2(X,Y,Z,90);
end

X0=XX(1,:);
Y0=XX(2,:);
Z0=XX(3,:);
