% Modified on 03/06/2023 fix the following:
%   1. match size of "radar_heatmap" to [N_rho,N_phi,N_theta]
%   2. automatically get two center rows of data with 2ss mode for different number of snapshots 
%%%--------------------------------------------------------
%%% This script was written based on SuperRF's
% https://bitbucket.org/embedded_intelligence/superrf_dataset/src/master/SuperRF_Data/
%%% -------------------------------------------------------
%%% The script will do the following jobs:
%   1-1. read raw adc data from the binary files, organize them into a 3d array
%   1-2. generate intensity/energy map from data in 1, in spherical coord.
%   2. convert the map from spherical coord. to Cartesian coord.
%%% -------------------------------------------------------
%%% Useful mmWave documentation:
%       swra581b.pdf -- Mmwave Radar Device ADC Raw Data Capture
%       spruij4a.pdf -- DCA1000EVM Data Capture Card
%%% -------------------------------------------------------
close all; clear; clc; 
format shortg;
%% preparation
Fs = 2047e3; % sampling frequency
sweep_slope = 29.982e+12; % Hz/s
FOV_phi = [30,150]; FOV_theta = [75,105]; 

full_or_2ss = "2ss";%"full"; % export mode: "full"(full-scale) or "2ss"(2-snapshot)
numFrames = 1;
numADCSamples = 256;
numChirp = 32; % number of chirps per frame 
numSnap = 64; % #rows, #snapshots, #vertical scans
numHz = 3; % #horizontal scans, 3*8=24 antennas for ISK
numTX = 2; numRX = 4;

rootaddr = ".\"; 
loadaddr = strcat(rootaddr); % sar data root folder
saveaddr = strcat(rootaddr,"heatmap\");

% sensor is at origin(0,0,0)
scene_lim = [-1, 1; 1, 3; -1, 1]; % scene boundary [x1,x2; y1,y2; z1,z2] (meters)
N_phi = 64; N_rho = 256; N_theta = 64; % dimension of spherical coord.
N_x= 64; N_y = 256; N_z = 64; % dimension of Cartesian coord.
threshold_factor = 0.2; % (<= 20% is suggested) for filtering
% values below (max intensity value*thershold_factor) will be dropped

% ------------------ modify parameters above as you need --------------- %

if ~exist(saveaddr, 'dir')
    mkdir(saveaddr);
    disp(strcat("Created directory: '", saveaddr,"'"));
end
%rf_complex = zeros(numSnap,numTX*numRX*numHz,numADCSamples); % 64x24x256
rf_complex_h1 = zeros(numSnap,numTX*numRX,numADCSamples); % 64x8x256
rf_complex_h2 = zeros(numSnap,numTX*numRX,numADCSamples);
rf_complex_h3 = zeros(numSnap,numTX*numRX,numADCSamples); 

%%% ================================================================== %%%
disp(" <--- Start ---> ") 
clk = clock; disp(clk);

%% read raw data and generate 3d heatmaps
% if you have data for more than one view, add an outer loop for it
for fm = 1:numFrames
    disp(strcat("frame # ",num2str(fm)))
    f = waitbar(0,'Loading data...');
    for hz = 1:numHz
        for sn = 1:numSnap
   %%% ------------------------- read data ------------------------- %%%
            waitbar(((hz-1)*numSnap+sn)/(numHz*numSnap),f,strcat("Loading hz-sn: ", num2str(hz),"-",num2str(sn)));
            %disp(strcat("hz-sn: ", num2str(hz),"-",num2str(sn)))
            addr = strcat(loadaddr,"\",num2str(hz),"\",num2str(sn),"\adc_data.bin");
            originData = readDCA1000(addr, numADCSamples); 

            numST = numADCSamples*numChirp; % total #samples for one TX
            framedata = originData(:,numST*numTX*(fm-1)+1 : numST*numTX*fm); % data in one frame

            TX1 = framedata(:,numST*0+1:numST*1);
            TX3 = framedata(:,numST*1+1:numST*2);

            % first chirp -> 256x1, turn them into column vectors
            T1R1 = reshape(TX1(1,1:numADCSamples),[numADCSamples,1]); 
            T1R2 = reshape(TX1(2,1:numADCSamples),[numADCSamples,1]);
            T1R3 = reshape(TX1(3,1:numADCSamples),[numADCSamples,1]);
            T1R4 = reshape(TX1(4,1:numADCSamples),[numADCSamples,1]);      
            T3R1 = reshape(TX3(1,1:numADCSamples),[numADCSamples,1]);
            T3R2 = reshape(TX3(2,1:numADCSamples),[numADCSamples,1]);
            T3R3 = reshape(TX3(3,1:numADCSamples),[numADCSamples,1]);
            T3R4 = reshape(TX3(4,1:numADCSamples),[numADCSamples,1]);

            % range-fft(sin(phase angle(rawdata))), 256x1
            T1R1r = fft(sin(angle(T1R1)));
            T1R2r = fft(sin(angle(T1R2)));
            T1R3r = fft(sin(angle(T1R3)));
            T1R4r = fft(sin(angle(T1R4)));
            T3R1r = fft(sin(angle(T3R1)));
            T3R2r = fft(sin(angle(T3R2)));
            T3R3r = fft(sin(angle(T3R3)));
            T3R4r = fft(sin(angle(T3R4)));

            % chirps in one row/snapshot
            TRs = [T1R1r T1R2r T1R3r T1R4r T3R1r T3R2r T3R3r T3R4r].'; % 8x256
            switch hz
                case 1
                    rf_complex_h1(sn,:,:) = TRs;
                case 2
                    rf_complex_h2(sn,:,:) = TRs;
                case 3
                    rf_complex_h3(sn,:,:) = TRs;
            end
        end
    end
    % flip over elements in each column (1st snapshot is at the bottom)
    rf_complex_h1 = flip(rf_complex_h1,1); % 64x8x256
    rf_complex_h2 = flip(rf_complex_h2,1); 
    rf_complex_h3 = flip(rf_complex_h3,1); 

    % h1 is the left-most column when facing towards the front of radar
    rf_complex = [rf_complex_h1,rf_complex_h2,rf_complex_h3]; % 64x24x256
    waitbar(1,f,'Finish loading data'); pause(1);close(f);
    disp('done intiallizing')

 %%% ------------------ generate 3d intensity maps ------------------ %%%
    disp('generating 3d intensity maps in spherical')
    % Hawkeye's data structure (range*azimuth*elevation)
    radar_heatmap = zeros(N_rho,N_phi,N_theta); 
    while true
        if full_or_2ss == "full"
            disp("Export mode: full-scale");
            rf_data_full = rf_complex; % full scale
            break;
        elseif full_or_2ss == "2ss"
            disp("Export mode: 2-snapshot");
            rf_data_2ss = rf_complex(ceil(numSnap/2):ceil(numSnap/2)+1,:,:); % 2 snapshots
            break;
        else
            disp("Incorrect export mode!")
            full_or_2ss = input("Re-enter (full or 2ss): ",'s');
        end 
    end
    
    numHzAtn = size(rf_complex, 2); % 24 horizontal virtual antenna
    fft_data = zeros(numSnap,numSnap);

    for range_idx = 1:numADCSamples
        if full_or_2ss == "full"
            fft_data(:,1:numHzAtn) = rf_data_full(:,:,range_idx); % full-scale
        else
            fft_data(1:2,1:numHzAtn) = rf_data_2ss(:,:,range_idx);
        end
        all_fft = fft2(fft_data); % angle fft
        all_fft = fftshift(all_fft);
        all_fft = flip(flip(all_fft,2));
        temp = abs(all_fft); 
        temp2 = temp.'; %(elv,azi) -> (azi,elv)
        radar_heatmap(range_idx,:,:) = temp2;
    end

    save(strcat(saveaddr,"radar_heatmap_sph_",num2str(fm),".mat"), 'radar_heatmap');
    disp(strcat('finished frame: ',num2str(fm))) 
    clk = clock; disp(clk);
end

heatmap_sph = radar_heatmap;
%heatmap_sph_scaled = scale_heatmap_color(heatmap_sph);
show_heatmap2d_sph(heatmap_sph);


%%% ================================================================== %%%


%% Convert spherical coordinates to Cartesian (within a specific space)

% convert points into cartesian coordinates
[x_ct,y_ct,z_ct] = sph2cart_pts(N_phi,N_rho,N_theta,Fs,sweep_slope,FOV_phi,FOV_theta);
ct_coord = [x_ct,y_ct,z_ct];

% match intensity values to corresponding spherical voxel center
radar_heat = matchHeat(heatmap_sph,N_phi,N_rho,N_theta);

% filtering: get points whose intensities are larger than the threshold
threshold = max(max(max(heatmap_sph)))*threshold_factor;

idx_heat_fliter = find(radar_heat >= threshold);
points_selected = zeros(size(idx_heat_fliter,1),3);
heat_selected = zeros(size(idx_heat_fliter,1),1);
for i = 1:size(idx_heat_fliter,1)
    points_selected(i,:) = ct_coord(idx_heat_fliter(i),:);
    heat_selected(i) = radar_heat(idx_heat_fliter(i));
end

% converting 3d heatmap into Cartesian
% assign value to the nearest Cartesian voxel directly from spherical voxels
disp('generating 3d intensity maps in Cartesian')
heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,points_selected,heat_selected);

save(strcat(saveaddr,"radar_heatmap_cart_",num2str(fm),".mat"), 'heatmap_ct');
disp("conversion: done");
clk = clock; disp(clk);


heatmap_ct_scaled = scale_heatmap_color(heatmap_ct);
show_heatmap2d_cart(heatmap_ct_scaled, scene_lim);
show_slice_heat_3d(heatmap_ct_scaled,scene_lim);

%{
% logarithmic scale
heatmap_ct = heatmap_ct + 1.0;
heatmap_ct = log(heatmap_ct);

% post filtering
maxheat = max(max(max(heatmap_ct)));
filter_idx = find(heatmap_ct < maxheat*0.675);
heatmap_ct(filter_idx) = 0;  
%}

%%% ================================================================== %%%


%%







                %%% --------------------------------- %%%
                %%% ---------- LEAVE BLANK ---------- %%%
                %%% --------------------------------- %%%

                
                

                
                
                
                
                
%% functions
% read raw adc data from binary files 
function [retVal] = readDCA1000(fileName, numADCSamples)
    %% global variables
    % change based on sensor config
    %numADCSamples = 256; % number of ADC samples per chirp
    numADCBits = 16; % number of ADC bits per sample
    numRX = 4; % number of receivers
    numLanes = 2; % do not change. number of lanes is always 2
    isReal = 0; % set to 1 if real only data, 0 if complex data0
    %% read file
    % read .bin file
    fid = fopen(fileName,'r');
    adcData = fread(fid, 'int16');
    % if 12 or 14 bits ADC per sample compensate for sign extension
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
    end
    fclose(fid);
    fileSize = size(adcData, 1);
    %disp(fileSize)
    %disp(size(adcData))
    
    % real data reshape, filesize = numADCSamples*numChirps
    if isReal
        numChirps = fileSize/numADCSamples/numRX;
        LVDS = zeros(1, fileSize);
        %create column for each chirp
        LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
        %each row is data from one chirp
        LVDS = LVDS.';
    else
        % for complex data
        % filesize = 2 * numADCSamples*numChirps
        numChirps = fileSize/2/numADCSamples/numRX;
        %disp("numChirps ="+numChirps)

        LVDS = zeros(1, fileSize/2);
        % Each LVDS lane contains 2 bytes, either 2I or 2Q. 
        %combine real and imaginary part into complex data
        %read in file: 2I is followed by 2Q
        counter = 1;
        for i=1:4:fileSize-1
            LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
            LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
            counter = counter + 2;
        end
        
        %disp("LVDS size ="+size(LVDS));
        
        % create column for each chirp
        LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
        %disp("LVDS="+size(LVDS));
        %each row is data from one chirp
        LVDS = LVDS.';
        %disp("LVDS="+size(LVDS));
    end
    %organize data per RX
    adcData = zeros(numRX,numChirps*numADCSamples);
    for row = 1:numRX
        for i = 1: numChirps
            adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = ...
                LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
        end
    end
    % return receiver data
    %disp("adcData size =" + size(adcData))
    retVal = adcData;
end

% convert points into cartesian coordinates
function [x_ct,y_ct,z_ct] = sph2cart_pts(N_phi,N_rho,N_theta,Fs,sweep_slope,FOV_phi,FOV_theta)
    c = 3e8; % speed of light 
    %BW = 3987.61e6; % Bandwidth = 3987.61MHz
    %Fs = 2047e3;
    %sweep_slope = 29.982e+12; % Hz/s
    %rho_res = c/2/BW; % in m, range resolution, 3.8cm for 3987.61MHz
    rho_min = 0; rho_max = Fs*c/2/sweep_slope; % range, m
    %phi_min = 30; phi_max = 150; % azimuth, deg
    %theta_min = 75; theta_max = 105; % elevation, deg
    phi_min = FOV_phi(1); phi_max = FOV_phi(2); % azimuth, deg
    theta_min = FOV_theta(1); theta_max = FOV_theta(2); % elevation, deg
    % get corresponding positions of each cell
    ticks_phi = linspace(phi_min,phi_max,N_phi); % azimuth axis of the output radar heatmap in degree
    ticks_theta = linspace(theta_min,theta_max,N_theta);
    ticks_rho = linspace(rho_min,rho_max,N_rho);
    ticks_phi = ticks_phi/180*pi;
    ticks_theta = ticks_theta/180*pi;
    % to use sph2cart, theta is the elevation angle from x-y plane
    ticks_theta_top = pi/2 - ticks_theta(1:32);
    ticks_theta_bottom = -(ticks_theta(33:64) - pi/2);
    ticks_theta = [ticks_theta_top,ticks_theta_bottom];

    N_cell = N_phi*N_theta*N_rho;
    sp_coord = zeros(N_cell,3);
    idx_sp_coord = 1;
    for idx_rho = 1:N_rho
        for idx_phi = 1:N_phi
            for idx_theta = 1:N_theta
                sp_coord(idx_sp_coord,:) = [ticks_rho(idx_rho),ticks_phi(idx_phi),ticks_theta(idx_theta)];
                idx_sp_coord = idx_sp_coord + 1;
            end
        end
    end
    [x_ct,y_ct,z_ct] = sph2cart(sp_coord(:,2),sp_coord(:,3),sp_coord(:,1)); % [x,y,z] = sph2cart(azimuth,elevation,r)
end

% match intensity values to corresponding points
function radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta)
    N_cell = N_phi*N_theta*N_rho;
    radar_heat = zeros(N_cell,1);
    idx_sp_coord = 1;
    for idx_rho = 1:N_rho
        for idx_phi = 1:N_phi
            for idx_theta = 1:N_theta
                radar_heat(idx_sp_coord) = heatmap(idx_rho,idx_phi,idx_theta); 
                idx_sp_coord = idx_sp_coord + 1;
            end
        end
    end
end

% match the heatmap into cartesian coordinates
function heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,pts,radar_heat)
    x_min = scene_lim(1,1); x_max = scene_lim(1,2);
    y_min = scene_lim(2,1); y_max = scene_lim(2,2);
    z_min = scene_lim(3,1); z_max = scene_lim(3,2);
    xx = linspace(x_min,x_max,N_x);
    yy = linspace(y_min,y_max,N_y);
    zz = linspace(z_min,z_max,N_z);
    % create a meshgrid and assign points to corresponding cubes
    [X,Y,Z] = meshgrid(xx,yy,zz);
    grid_centers = [X(:),Y(:),Z(:)];
    %clss = knnsearch(grid_centers,[x_ct,y_ct,z_ct]); % classification
    % pts = [x_ct,y_ct,z_ct];
    clss = knnsearch(grid_centers,pts,'K',1); % classification
    local_stat = @(x)mean(x); % defintion of local statistic
    %class_stat = accumarray(clss,radar_heat,[numr*numc*256 1],local_stat); % data_grouping
    class_stat = accumarray(clss,radar_heat,[N_x*N_y*N_z 1],local_stat); % data_grouping
    heatmap_ct  = reshape(class_stat , size(X)); % 3D reshaping
    
    for id_y = 1:N_y
        for id_x = 1:N_x
            for id_z = 1:N_z
                if heatmap_ct(id_y,id_x,id_z) == 0
                    if id_x == 1
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/2;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/2;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        end
                    elseif id_x == N_x
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/2;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/2;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        end
                    else
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/3;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/4;
                        end
                    end
                end
            end
        end
    end
    
end

% plot 2d radar heatmaps in spherical coordinates
function f = show_heatmap2d_sph(heatmap)
    f = figure(); sgtitle("Spherical 3-view drawing"); font_size = 8;
    az_tick = linspace(1,64,9); rg_tick = linspace(1,256,11); el_tick = linspace(1,64,7);
    az_label = {'-60','-45','-30','-15','0','15','30','45','60'};
    rg_label = {'0','1','2','3','4','5','6','7','8','9','10'};
    el_label = {'-15','-10','-5','0','5','10','15'};
    
    % Visulize the radar heatmap top view
    radar_heatmap_top = squeeze(max(heatmap,[],3));
    subplot(131); imagesc(radar_heatmap_top);  title("Top"); 
    set(gca,'XDir','reverse'); set(gca,'YDir','normal');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Azimuth (deg)'); ylabel('Range (m)'); set(gca,'FontSize',font_size);
    xticks(az_tick); xticklabels(az_label)
    yticks(rg_tick); yticklabels(rg_label)
    xlim([1,64]),ylim([1,256])
    
    % Visulize the radar heatmap front view
    radar_heatmap_front = squeeze(max(heatmap,[],1));
    subplot(132); imagesc(radar_heatmap_front.');  title("Front");
    set(gca,'XDir','reverse');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)'); set(gca,'FontSize',font_size);
    xticks(az_tick); xticklabels(az_label)
    yticks(el_tick); yticklabels(el_label)
    xlim([1,64]),ylim([1,64])
    
    % Visulize the radar heatmap side view
    radar_heatmap_side = squeeze(max(heatmap,[],2));
    subplot(133); imagesc(radar_heatmap_side.'); title("Side"); 
    %set(gca,'XDir','reverse');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Range (m)'); ylabel('Elevation (deg)'); set(gca,'FontSize',font_size);
    xticks(rg_tick); xticklabels(rg_label)
    yticks(el_tick); yticklabels(el_label)
    xlim([1,256]),ylim([1,64])
    
    fp = f.Position;
    f.Position = [fp(1),fp(2),fp(3)*1.5,fp(3)/2];
end

% plot 2d radar heatmaps in cartesian coordinates
function f = show_heatmap2d_cart(heatmap, scene_lim)
    f = figure(); sgtitle("Cartesian 3-view drawing"); font_size = 8;
    x_tick = linspace(1,64,11); y_tick = linspace(1,256,11); z_tick = linspace(1,64,11);
    x_label = num2cell(linspace(scene_lim(1,1),scene_lim(1,2),11));%{'-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'};
    y_label = num2cell(linspace(scene_lim(2,1),scene_lim(2,2),11));%{'0','1','2','3','4','5','6','7','8','9','10'};
    z_label = num2cell(linspace(scene_lim(3,1),scene_lim(3,2),11));%{'-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'};
    
    % Visulize the radar heatmap top view
    radar_heatmap_top = squeeze(max(heatmap,[],3));
    subplot(131); imagesc(radar_heatmap_top); title("Top");
    %set(gca,'XDir','reverse');
    set(gca,'XDir','normal');set(gca,'YDir','normal');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('x (m)'); ylabel('y (m)'); set(gca,'FontSize',font_size);
    xticks(x_tick); xticklabels(x_label)
    yticks(y_tick); yticklabels(y_label)
    xlim([1,64]),ylim([1,256])
    
    % Visulize the radar heatmap front view
    radar_heatmap_front = squeeze(max(heatmap,[],1));
    subplot(132); imagesc(radar_heatmap_front.'); title("Front");
    set(gca,'XDir','normal'); set(gca,'YDir','normal');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('x (m)'); ylabel('z (m)'); set(gca,'FontSize',font_size);
    xticks(x_tick); xticklabels(x_label)
    yticks(z_tick); yticklabels(z_label)
    xlim([1,64]),ylim([1,64])
    
    % Visulize the radar heatmap side view
    radar_heatmap_side = squeeze(max(heatmap,[],2));
    subplot(133); imagesc(radar_heatmap_side.'); title("Side");
    set(gca,'XDir','normal');set(gca,'YDir','normal');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('y (m)'); ylabel('z (m)'); set(gca,'FontSize',font_size);
    xticks(y_tick); xticklabels(y_label)
    yticks(z_tick); yticklabels(z_label)
    xlim([1,256]),ylim([1,64])  
    
    fp = f.Position;
    f.Position = [fp(1),fp(2),fp(3)*1.5,fp(3)/2];
end

% plot 3d heatmap slice in cartesian
function f = show_slice_heat_3d(heatmap_ct,scene_lim)
    m = size(heatmap_ct, 2); % az
    n = size(heatmap_ct, 3); % el
    l = size(heatmap_ct, 1); % rg
    yi = linspace(1,l,l); % range
    xi = linspace(1,m,m); % azimuth
    zi = linspace(1,n,n); % elevation
    [XX,YY,ZZ] = meshgrid(xi,yi,zi); % [l * m * n]

    xslice = 1:m;    % location of y-z planes
    yslice = 1:l;    % location of x-z plane
    zslice = 1:n;    % location of x-y planes
    
    f = figure('visible', 'on');
    h = slice(XX,YY,ZZ,heatmap_ct,xslice,yslice,zslice);
    title("Cartesian 3D intensity map")
    xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
    xlim([0 64]); ylim([0 256]); zlim([0,64]);
    ticklables1 = num2cell(linspace(scene_lim(1,1),scene_lim(1,2),11));
    ticklables2 = num2cell(linspace(scene_lim(2,1),scene_lim(2,2),11));
    ticklables3 = num2cell(linspace(scene_lim(3,1),scene_lim(3,2),11));

    xt = linspace(1,64,11); xticks(xt); 
    xticklabels(ticklables1);
    yt = linspace(0,256,11); yticks(yt); 
    yticklabels(ticklables2);
    zt = linspace(0,64,11); zticks(zt); 
    zticklabels(ticklables3);
    %set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');
    set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','flat');
    %set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','0.01');
    % set transparency to correlate to the data values.
    alpha('color'); %alpha(h, 0.01);
    colorbar; colormap jet;
end

% scale down the intensity for better plots
function heatmap_scaled = scale_heatmap_color(heatmap)
    div = floor(log10(max(max(max(heatmap))))) - 2;
    heatmap_scaled = heatmap/(10^div);
    heatmap_scaled = log(heatmap_scaled+1);

    threshold = max(max(max(heatmap_scaled)))*0.01;
    heatmap_scaled(heatmap_scaled<=threshold)=0;
end