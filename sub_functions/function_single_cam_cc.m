%Function four_cams2 runs collinearity equations for 4-camera setup

%Function four_cams2 takes as its arguments:
   
    % - x (an 8-by-n matrix containing pixel coordinates of the n points in the four cameras in the order
        %[x1,y1,x2,y2,x3,y3,x4,y4])
    % - the 3-by-18+n matrix chi, concatenating:
        % - Rodrigues vectors for cameras 2-4 (:,1:3)
        % - Pose vectors for cameras 2-4 (:,4-6)
        % - A 3-by-3 matrix for each of cameras 1-4 containing:
            % - the focal length in horizontal (1,1) and vertical (2,1)
            % pixels
            % - the principal point coordinates (3,1) to (1,2)
            % - the tangential distortion coefficients P1 and P2 (2:3,2)
            % - the radial distortion coefficents K1, K2, K3 (:,3)    
    % - X (a 3-by-n matrix containing the X, Y, Z coordinates of the n
        %points in the axis system of camera 1)
   
            
function F=function_single_cam_cc(xx,chi,XX,kc)
% 
% save function_sing_cam_no_cc_temp
% return
% clear
% close all
% load function_sing_cam_no_cc_temp
% return

om_grids=chi(:,2:size(XX,2)+1);
T_grids=chi(:,size(XX,2)+2:end);

X=[];

for i=1:size(XX,2)
    X_temp=XX{i};
    X=[X rodrigues(om_grids(:,i))*X_temp+T_grids(:,i)*ones(1,size(X_temp,2))];
end

x_all=xx(1:2:end,:);
y_all=xx(2:2:end,:);

%Number of points
n=size(xx,2);

cc=chi(2:3,1);

dx=x_all-cc(1);
dy=y_all-cc(2);

r=dx.^2+dy.^2;

cx=chi(1,1);
cy=chi(1,1);

R=rodrigues(zeros(3,1));
T=zeros(3,n);

K=kc(1:3);
P=kc(4:5);
 
F=[cx*(R(1,:)*(T-X))./(R(3,:)*(T-X))-dx.*(ones(1,n)+K(1)*r+K(2)*r.^2+K(3)*r.^3)+P(1)*(r+2*dx.^2)+2*P(2)*dx.*dy;
   cy*(R(2,:)*(T-X))./(R(3,:)*(T-X))-dy.*(ones(1,n)+K(1)*r+K(2)*r.^2+K(3)*r.^3)+P(2)*(r+2*dy.^2)+2*P(1)*dx.*dy];

Fnan=isnan(F);
F(Fnan)=0;