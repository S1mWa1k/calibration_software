function [x_cam,XX_cam,dx,dy,nx0,ny0]=findCheckerboards(camera,base_folder,date,folder,draw_images,boardSize,I_check,imagePoints,ctrl_marker_pos)

% close all
% clear
% 
% base_folder='D:\Cunyi_bumblebees';
% date='2025_06_19';
% folder='Calib_mraw';
% draw_images=1;
% 
% d=dir(fullfile(base_folder,date,folder,'*.cihx'));
% 
% Cameras={d.name};

% Cameras={'Camera_01';
%          'Camera_02';
%          'Camera_03';
%          'Camera_06'};

% for c=1:length(Cameras)

    % clearvars -except c Cameras base_folder date folder draw_images
    % close all
% 
    % camera=Cameras{c};

    disp(['Calibrating ' camera])
    
    % draw_images=1;
    save_images=0;
    
    if save_images==1 && ~exist([folder ' tifs'],'dir')
        mkdir([folder ' tifs'])
    end
    
    video=quick_parse_cihx(fullfile(base_folder,date,folder,camera));
    n_im=video.TotalFrame;
    
    I=zeros(video.ImageHeight,video.ImageWidth,1,n_im,'uint16');
    warning off
    
    boardSize_real=boardSize-1;
    % boardSize_real=[8 9]; % size of board in rows x columns
    % ctrl_marker_pos=[4 4]; % position of white ctrl marker (rows x column)
        
    px0=reshape(imagePoints(:,1),boardSize(1)-1,[]);
    py0=reshape(imagePoints(:,2),boardSize(1)-1,[]);
        
    % tic
    % imagesc2(I_check); hold on
    % plot(px0,py0,'.-','markersize',15,'linewidth',2)
    
    % for i=1:numel(px0)
    % 
    %     text(px0(i)+10,py0(i)-20,num2str(i),'Color','r')
    % 
    % end
    % toc
    % return
    
    for f=1:n_im
    
        x_cam0=NaN(boardSize_real);
        y_cam0=NaN(boardSize_real);
    
        x_cam{f}=[x_cam0(:) y_cam0(:)]';
    
        im=readmraw(video,f);
        % im=imread('new_checkerboard.tif');
        
         % bw = imbinarize(im);
        % im2=imgaussfilt(im,2);
    
        [imagePoints,boardSize] = detectCheckerboardPoints(im,'MinCornerMetric',0.2);
    
        if min(boardSize)>5 && max(boardSize)<=max(boardSize_real+1)
            imagePoints(imagePoints(:,1)<10,:)=NaN;
            imagePoints(imagePoints(:,2)<10,:)=NaN;
            imagePoints(imagePoints(:,1)>size(im,2)-9,:)=NaN;
            imagePoints(imagePoints(:,2)>size(im,1)-9,:)=NaN;
    
            % imagePoints(imagePoints(:,1)<100,:)=NaN;
            % imagePoints(imagePoints(:,2)<100,:)=NaN;
            % imagePoints(imagePoints(:,1)>size(im,2)-99,:)=NaN;
            % imagePoints(imagePoints(:,2)>size(im,1)-99,:)=NaN;
            
            % imagesc2(im); hold on
            % plot(imagePoints(:,1),imagePoints(:,2),'.g','markersize',15);
            
            if sum(~isnan(imagePoints(:,1)))>20
                
                px=reshape(imagePoints(:,1),boardSize(1)-1,[]);
                py=reshape(imagePoints(:,2),boardSize(1)-1,[]);
                
                % get mid-points of each square
                px_mid=0.5*(px(:,1:end-1)+px(:,2:end));
                py_mid=0.5*(py(:,1:end-1)+py(:,2:end));
            
                px_mid2=0.5*(px_mid(1:end-1,:)+px_mid(2:end,:));
                py_mid2=0.5*(py_mid(1:end-1,:)+py_mid(2:end,:));
            
                px_mid2(isnan(py_mid2))=NaN;
                py_mid2(isnan(px_mid2))=NaN;
        
                % imagesc2(im); hold on
                % plot(px,py,'.-','markersize',15,'linewidth',2)
                % plot(px_mid2,py_mid2,'.r','markersize',15)
                % return
        
                % blur image
                imb=imgaussfilt(im,5);
            
                % get values of square centres
                px_mid22=px_mid2;
                py_mid22=py_mid2;
                
                px_mid22(isnan(px_mid2))=1;
                py_mid22(isnan(py_mid2))=1;
        
                pind_mid2=sub2ind(size(im),round(py_mid22),round(px_mid22));
                
                point_int=double(imb(pind_mid2));
                point_int(isnan(py_mid2))=NaN;
        
                % imagesc2(im); hold on
                % plot(px_mid2,py_mid2,'.r','markersize',15)
                % return
        
                point_int=point_int-mean(point_int(:),'omitnan');
            
                % create blank mesh
                BW01=repmat([-1 1;1 -1],ceil(0.5*size(point_int,1)),ceil(0.5*size(point_int,2)));
            
                BW01=BW01(1:size(point_int,1),1:size(point_int,2));
                
                R = corrcoef(point_int,BW01,'rows','complete');
        
                % disp(abs(R(2)))
                
                if abs(R(2))>0.5
                
                    if R(2)<0
                        BW01=BW01*-1;
                    end
                
                    point_match=BW01.*point_int;
                
                    ctrl_points=find(point_match<0);
        
                    if length(ctrl_points)==2
                
                        if point_int(ctrl_points(1))<point_int(ctrl_points(2))
                    
                            ctrl_points=flip(ctrl_points);
                        end
                    
                        [py_ctrl_ind, px_ctrl_ind]=ind2sub(size(point_match),ctrl_points);
                        
                        px_ctrl=(px_mid2(ctrl_points));
                        py_ctrl=(py_mid2(ctrl_points));
            
                        %% now map found pounts onto original grid
                        rotated=0;
            
                        %% rotate grid lines if running wrong way
                        if diff(py_ctrl_ind)~=0
            
                            px=px';
                            py=py';
            
                            % px=flip(px,2);
                            % py=flip(py,2);
                            
                            temp_offset=py_ctrl_ind;
                            py_ctrl_ind=px_ctrl_ind;
                            px_ctrl_ind=temp_offset;
                            
                            rotated=1;
                        end
                        
                        if ~any(size(px)-boardSize_real>0) && ~any(px_ctrl_ind>5) && ~any(size(px,2)-px_ctrl_ind>5) && ~any(py_ctrl_ind>4) && ~any(size(py,1)-py_ctrl_ind>4)
    
                            % imagesc2(im); hold on
                            % plot(px,py,'.-','markersize',15,'linewidth',2)
                            % % plot(px_mid2,py_mid2,'.r','markersize',15)
                            % plot(px_ctrl(1),py_ctrl(1),'.r','markersize',16);
                            % plot(px_ctrl(2),py_ctrl(2),'.g','markersize',16);
                            % return
                            
                            if diff(py_ctrl_ind)==0
                
                                if diff(px_ctrl_ind)==1
                                    py_ctrl_offset=ctrl_marker_pos(1)-py_ctrl_ind(1)+1;
                                    px_ctrl_offset=ctrl_marker_pos(2)-px_ctrl_ind(1)+1;
                                    
                                    x_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1,px_ctrl_offset:px_ctrl_offset+size(px,2)-1)=px;
                                    y_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1,px_ctrl_offset:px_ctrl_offset+size(px,2)-1)=py;
                                
                                    flipped=0;
                                    
                                    % if mean(py(1,:)-py(end,:),'omitnan')>0
                                    %     x_cam0=flip(x_cam0,1);
                                    %     y_cam0=flip(y_cam0,1);
                                    % end
                                    
                                else
                
                                    py_ctrl_offset=ctrl_marker_pos(1)-(size(px,1)-py_ctrl_ind(1))+1;
                                    px_ctrl_offset=ctrl_marker_pos(2)-(size(px,2)-px_ctrl_ind(1))+1;
                
                                    % rotate 180 deg
                                    px=flip(px,2);
                                    py=flip(py,2);
                
                                    px=flip(px,1);
                                    py=flip(py,1);
                                                        
                                    x_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1,px_ctrl_offset:px_ctrl_offset+size(px,2)-1)=px;
                                    y_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1,px_ctrl_offset:px_ctrl_offset+size(px,2)-1)=py;
                
                                    
                                    % if mean(py(1,:)-py(end,:),'omitnan')>0
                                    %     x_cam0=flip(x_cam0,1);
                                    %     y_cam0=flip(y_cam0,1);
                                    % end
                
                                    % x_cam0=flip(x_cam0,2);
                                    % y_cam0=flip(y_cam0,2);
                                    % 
                                    % x_cam0=flip(x_cam0,1);
                                    % y_cam0=flip(y_cam0,1);
                
                                    flipped=1;
                                end
                            else
                
                                moo=1
                                return
                            end
                            
                            % ensure correct orientation
                            P0 = [px_ctrl(2), py_ctrl(2)];
                            P1 = [px_ctrl(1), py_ctrl(1)];
                            P2 = [px(:,1), py(:,1)];
            
                            angle4=NaN(size(px,1),1);
            
                            for i=1:size(P2,1)
                                angle4(i) = angle3Points(P1, P0, P2(i,:));
                                
                                % n1 = (P2(i,:) - P0) / norm(P2(i,:) - P0);  % Normalized vectors
                                % n2 = (P1 - P0) / norm(P1 - P0);
                                % angle3(i) = atan2(norm(det([n2; n1])), dot(n1, n2));
                            end
            
                            angle4=unwrap(angle4);
            
                            if mean(diff(angle4),'omitnan')>0
            
                                x_cam0=flip(x_cam0,1);
                                y_cam0=flip(y_cam0,1);
            
                            end
                            % return
        
                            x_cam{f}=[x_cam0(:) y_cam0(:)]';
        
                            if draw_images==1
            
                                subplot(1,2,1)
                                hold off
                                imagesc2(im); hold on
                                % plot(px,py,'.-','markersize',15,'linewidth',2)
                                plot(x_cam0,y_cam0,'.-','markersize',15,'linewidth',2)
                                % plot(px_mid2,py_mid2,'.r','markersize',15)
                                plot(px_ctrl(1),py_ctrl(1),'.r','markersize',16);
                                plot(px_ctrl(2),py_ctrl(2),'.g','markersize',16);
                    
                                for i=1:numel(x_cam0)
                    
                                    text(x_cam0(i)+10,y_cam0(i)-20,num2str(i),'Color','r')
                    
                                end
                    
                                title([num2str(f) ' / ' num2str(n_im)])
                                axis off
                    
                                subplot(1,2,2)
                                hold off
                                px00=px0;
                                py00=py0;
                    
                                px00(isnan(x_cam0))=NaN;
                                py00(isnan(x_cam0))=NaN;
                                
                                imagesc2(I_check); hold on
                                plot(px00,py00,'.-','markersize',15,'linewidth',2)
                                axis off

                                title(['Camera ' camera])
                    
                                set(gcf,'color','w','position',[100 100 900 450])
                        
                                for i=1:numel(px00)
                    
                                    text(px00(i)+10,py00(i)-20,num2str(i),'Color','r')
                    
                                end
                                
                                drawnow
                                % return
                
                                if save_images==1
                             
                                    set(gcf,'PaperPositionMode', 'auto')
                                    print('-dpng','-r0',fullfile([folder ' tifs'],[folder ' f' num2str(f) '.png']));
                                end
        
                            end
                        else
    
                        end
        
                    else
        
        
                    end
                
                else
                
                end
            
            else
    
            end
    
        else
        
        end
    
    end
    
    ny0=boardSize_real(2);
    nx0=boardSize_real(1);
    nd=16;
    dx=video.ImageWidth;
    dy=video.ImageHeight;
    
    XX_cam=[reshape(repmat(0:ny0-1,nx0,1),1,[]);repmat(0:nx0-1,1,ny0);ones(1,nx0*ny0)];
    XX_cam(1:2,:)=XX_cam(1:2,:)*nd;
    XX_cam=repmat({XX_cam},1,max(n_im));
    
    % save(fullfile(base_folder,date,[folder '_processed' '_' camera]),'x_cam','XX_cam','nx0','ny0','nd','dx','dy','n_im')
   close
% end