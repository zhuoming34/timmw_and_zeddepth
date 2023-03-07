%%% This script scales real depth images to a size with specified boundaries
close all; clear; clc;
format shortg;

rootaddr = ".\";
loadaddr = strcat(rootaddr);
saveaddr = strcat(rootaddr,"processed images\");

numImgFrame = 1; % number of frames to be exported
% the space selected for scaling depth images
% 4 corner points in (x,y,z), [lower-left, upper-left, l-right, u-right]
% y cannot be 0 !!! x_left = -x_right (symmetric)
scene_corners = [-1,1,-1; -1,1,1;  1,1,-1; 1,1,1]; 

% ------------------ Modify Parameters Above As You Need --------------- %
% ---------------------------------------------------------------------- %


%% Preparation
% basic parameters with 720P resolution from ZED
numc = 1280; numr = 720; 
rng_min = 0.1*1000; rng_max = 15.0*1000; % mm, min/max range of ZEDmini
focal_px = 700; % focal length in number of pixels
pxSize = 0.004; % pixel size, mm
focal_L = focal_px*pxSize/1000; % focal length 2.8mm -> 0.0028m
% ---------------------------------------------------------------------- %
map1 = jet; map2 = gray;
% ---------------------------------------------------------------------- %
% projection plane 
ppH = numc*pxSize/1000; ppV = numr*pxSize/1000; % size of plane, in m
px_l = -ppH/2; px_r = ppH/2; pz_b = -ppV/2; pz_t = ppV/2; % boundary: left, right, bottom, top

% scene boundary and their projections
sn_x = scene_corners(:,1); sn_y = scene_corners(:,2); sn_z = scene_corners(:,3);
sn_py = linspace(focal_L,focal_L,size(sn_y,1))'; % focal length or y values of projefction plane
sn_ratio = sn_y./sn_py; 
sn_px = -sn_x./sn_ratio;
sn_pz = -sn_z./sn_ratio;

% file management
subfolder1 = "1280x720\"; subfolder2 = "128x128\";
subfolder3 = "scaled\"; subfolder4 = "1280x720(colored)\";
checkFolders(saveaddr,subfolder1,subfolder2,subfolder3,subfolder4);

disp(" <--- Start ---> ") 
clk = clock; disp(clk);
%% Image processing
for img_idx = 1:numImgFrame
        if img_idx <= 9
            file_prefix = 'depth00000';
        else
            file_prefix = 'depth0000';
        end
        filename_mat = strcat(loadaddr,"mat files\",file_prefix,num2str(img_idx),'.mat');
        filename_img = strcat(loadaddr,"png files\",file_prefix,num2str(img_idx),'.png');
        depthmat = load(filename_mat).depimg;
        orgImg = imread(filename_img);

        depthmat(1,numc/2) = rng_min;
        depthmat(end,numc/2) = rng_max;
        depthmat(depthmat <= rng_min) = rng_min;
        depthmat(depthmat >= rng_max) = rng_max;                    

        % normalize pixel values to [0,1]
        I0 = depthmat;
        I = ( I0 - min(min(I0)) ) ./ ( max(max(I0)) - min(min(I0)) );
        
        % extend and crop to the same size of scene  
        [I2,NN_row,NN_col] = depScale(I,numr,numc,px_l,px_r,pz_t,pz_b,sn_px,sn_pz,pxSize);
        disp(strcat("Scaled size of img is: ",num2str(NN_row),"x",num2str(NN_col)));
        
        % saving images
        DepImg = abs(I - 1)*255; % original size 1280x720
        SclImg = abs(I2 - 1)*255; % extended
        SqrImg = imresize(SclImg, [128,128]); % squared extended

        % map1 = jet; map2 = gray;
        outputFileName1 = strcat(saveaddr,subfolder1,num2str(img_idx),'.png');
        imwrite(DepImg, map2, outputFileName1); 

        outputFileName2 = strcat(saveaddr,subfolder2,num2str(img_idx),'.png');
        imwrite(SqrImg, map2, outputFileName2); 
        
        outputFileName3 = strcat(saveaddr,subfolder3,num2str(img_idx),'.png');
        imwrite(SclImg, map1, outputFileName3); 
        
        outputFileName4 = strcat(saveaddr,subfolder4,num2str(img_idx),'.png');
        imwrite(DepImg, map1, outputFileName4);

        
        disp(strcat("Frame #",num2str(img_idx)," is done"));
        clk = clock; disp(clk);
end
disp('All finished')
clk = clock; disp(clk);

% plots
plot_depimg(DepImg,rng_min,rng_max,"1280x720");
plot_depimg(SclImg,rng_min,rng_max,"Scaled size");
plot_depimg(SqrImg,rng_min,rng_max,"128x128");
figure(); imshow(orgImg); title("origianlly from camera")

%%







                %%% --------------------------------- %%%
                %%% ---------- LEAVE BLANK ---------- %%%
                %%% --------------------------------- %%%

                
                

                
                
                
                
           
%% functions

function [I_scaled,NN_row,NN_col] = depScale(I,numr,numc,px_l,px_r,pz_t,pz_b,sn_px,sn_pz,pxSize)     
    % number of columns or rows to be cropped or extended
    N_col_r = ceil(abs((sn_px(3)-px_l))/(pxSize/1000));
    N_col_l = ceil(abs((sn_px(1)-px_r))/(pxSize/1000));
    N_row_b = ceil(abs((sn_pz(1)-pz_t))/(pxSize/1000));
    N_row_t = ceil(abs((sn_pz(2)-pz_b))/(pxSize/1000));

    % assume it is symmetric about x=0
    if sn_px(3) <= px_l || sn_px(1) >= px_r
        col_r = ones(numr,N_col_r); col_l = ones(numr,N_col_l); % extend
        Ic = [col_l, I, col_r];
        NN_col = numc + N_col_r + N_col_l; % new number of columns
    else
        Ic = I(:, N_col_l+1:end-N_col_r); % crop
        NN_col = numc - N_col_l - N_col_r;
    end

    if sn_pz(1) >= pz_t % bottom of img, positive on projection plane
        row_b = ones(N_row_b,NN_col); % extend
        Ib = [Ic;row_b];
        NN_row = numr + N_row_b; % new number of rows
    else
        Ib = Ic(1:end-N_row_b,:); % crop
        NN_row = numr - N_row_b;
    end

    if sn_pz(2) <= pz_b % top of img, negative on projection plane
        row_t = ones(N_row_t,NN_col); % extend
        It = [row_t;Ib];
        NN_row = NN_row + N_row_t;
    else
        It = Ib(N_row_t+1:end,:); % crop
        NN_row = NN_row - N_row_t;
    end
    
    I_scaled = It;
end

function checkFolders(saveaddr,subfolder1,subfolder2,subfolder3,subfolder4)
    if ~exist(saveaddr, 'dir')
        mkdir(saveaddr);
        disp(strcat("Created directory: '", saveaddr,"'"));
    end
    if ~exist(strcat(saveaddr,subfolder1), 'dir')
        mkdir(strcat(saveaddr,subfolder1));
        disp(strcat("Created directory: '", strcat(saveaddr,subfolder1),"'"));
    end
    if ~exist(strcat(saveaddr,subfolder2), 'dir')
        mkdir(strcat(saveaddr,subfolder2));
        disp(strcat("Created directory: '", strcat(saveaddr,subfolder2),"'"));
    end
    if ~exist(strcat(saveaddr,subfolder3), 'dir')
        mkdir(strcat(saveaddr,subfolder3));
        disp(strcat("Created directory: '", strcat(saveaddr,subfolder3),"'"));
    end
    if ~exist(strcat(saveaddr,subfolder4), 'dir')
        mkdir(strcat(saveaddr,subfolder4));
        disp(strcat("Created directory: '", strcat(saveaddr,subfolder4),"'"));
    end
end

function f = plot_depimg(img,rng_min,rng_max,myTitle)
        f = figure('visible','on');
        %imagesc(img); %axis equal;
        sgtitle(myTitle);
        image(img); axis image
        set(gca,'visible','off');
        colormap jet;
        c = colorbar; 
        %a = f.Position; 
        c.Label.String = 'range (m)';
        c.Ticks = linspace(1,256,16);
        c.TickLabels = num2cell(linspace(floor(rng_max/1000),floor(rng_min/1000),16));
        caxis([1,256]);
        c.Direction = 'reverse';
end