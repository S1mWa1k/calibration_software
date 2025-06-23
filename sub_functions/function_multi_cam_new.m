%Function function_multi_cam runs collinearity equations for 4-camera setup

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
   
            
function F=function_multi_cam_new(xx,chi,XX,fc,cc,kc,num_cam)
% % 
% save temp_chi_new
% return
% clear all
% load temp_chi_new

om_cams_2_ref=chi(:,1:num_cam);
T_cams_2_ref=chi(:,num_cam+1:2*num_cam);

om_grids_2_ref=chi(:,2*num_cam+1:2*num_cam+size(XX,2));
T_grids_2_ref=chi(:,2*num_cam+size(XX,2)+1:end);

%calculate 3xn matrix X from grid rotations and pose from chi
X=[];

for i=1:size(XX,2)
    X_temp=XX{i};
    X=[X rodrigues(om_grids_2_ref(:,i))*X_temp+T_grids_2_ref(:,i)*ones(1,size(X_temp,2))];
end

x_all=xx(1:2:end,:);
y_all=xx(2:2:end,:);

%Number of points
n=size(xx,2);

% dx_all=NaN(num_cam,n);
% dy_all=NaN(num_cam,n);
% 
% for i=1:num_cam
% 
%     dx_all(i,:)=x_all(i,:)-cc(1,i);
%     dy_all(i,:)=y_all(i,:)-cc(2,i);
% 
% end

for i=1:num_cam
    
    cx=fc(1,i);
    cy=fc(2,i);
    
    R=rodrigues(om_cams_2_ref(:,i))';
    T=repmat(T_cams_2_ref(:,i),1,n);
    
    K=kc(1:3,i);
    P=kc(4:5,i);

    dx=x_all(i,:)-cc(1,i);
    dy=y_all(i,:)-cc(2,i);
    
    % dx=dx_all(i,:);
    % dy=dy_all(i,:);

    dx_n=dx/cx;
    dy_n=dy/cy;

    r_n=sqrt(dx_n.^2+dy_n.^2);
        
    F(2*i-1:2*i,:)=[cx*(R(1,:)*(T-X))./(R(3,:)*(T-X))-(dx.*(ones(1,n)-(K(1).*r_n.^2+K(2).*r_n.^4+K(3).*r_n.^6)+2*P(1).*dx_n.*dy_n+P(2)*(r_n.^2+2*dx_n.^2)));
                    cy*(R(2,:)*(T-X))./(R(3,:)*(T-X))-(dy.*(ones(1,n)-(K(1).*r_n.^2+K(2).*r_n.^4+K(3).*r_n.^6)+2*P(2).*dx_n.*dy_n+P(1)*(r_n.^2+2*dy_n.^2)))];

end

Fnan=isnan(F);
F(Fnan)=0;