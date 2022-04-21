%% parameters for static features
pixel_resolution_static = .5; %degrees of visual angle
gauss_filter_size_static = 1; %degrees of visual angle, s.d. of Gaussian or radius of disk
filter_gauss = false; %true = gaussian, false = disk filter
OS_thres_frac = .5; %fraction contrast
contrast_thres_frac = .1;

%which features
features_static.local_contrast = true;
features_static.horizontal_orientation = true;
features_static.vertical_orientation = true;

%% parameters for motion features
pixel_resolution_motion = 0.2; %degrees
feature_spacing_motion = 2; %degrees
motion_feature_RF_center = 1; %degrees, could be fixed at spacing
motion_feature_RF_surround = 5; %degrees
motion_surround_weight = 2; %relative to center weight
OMS_threshold_frac = .1; %fraction contrast

%which features
features_motion.OMS = true;

%% Constants
video_degrees_visual_angle = 120;
contrast_max = 2^7; % for 8 bit image, positive or negative contrast.

%% Computed parameters
OMS_thres = round(OMS_threshold_frac * contrast_max);
contrast_thres = round(contrast_thres_frac * contrast_max);
OS_thres = round(OS_thres_frac * contrast_max);

%% Load video
video_fname = 'GameRecording_short.mov';
v = VideoReader(video_fname);
degrees_per_pixel = video_degrees_visual_angle / v.Width;

resize_factor_static = degrees_per_pixel / pixel_resolution_static;
h_static = round(v.Height*resize_factor_static);
w_static = round(v.Width*resize_factor_static);

resize_factor_motion = degrees_per_pixel / pixel_resolution_motion;
h_motion = round(v.Height*resize_factor_motion);
w_motion = round(v.Width*resize_factor_motion);

resize_factor_motion_feature = round(feature_spacing_motion / pixel_resolution_motion);
h_motion_feature = round(h_motion/resize_factor_motion_feature);
w_motion_feature = round(w_motion/resize_factor_motion_feature);

frames_static = zeros(h_static,w_static,v.NumFrames,'uint8');
frames_diff = zeros(h_motion,w_motion,v.NumFrames,'uint8');

i=1;
while hasFrame(v)
    curFrame = mean(readFrame(v), 3); %grayscale
    frames_static(:,:,i) = imresize(curFrame, resize_factor_static);
    if i>1        
        frames_diff(:,:,i) = imresize(curFrame - lastFrame, resize_factor_motion);
    end
    lastFrame = curFrame;
    i=i+1;
end

fprintf('Movie loaded and resampled (static) to: %d frames, %d x %d pixels, %d degrees per pixel\n', ...
    v.NumFrames, w_static, h_static, pixel_resolution_static);

fprintf('Movie loaded and resampled (motion) to: %d frames, %d x %d pixels, %d degrees per pixel\n', ...
    v.NumFrames, w_motion, h_motion, pixel_resolution_motion);

%TODO: load if these already exist
[h_static, w_static, frames] = size(frames_static);
[h_motion, w_motion, frames] = size(frames_diff);

%% subunits for motion computations
[posX, posY] = makeSubunitHexGrid([w_motion,h_motion], resize_factor_motion_feature);
N_motion_units = length(posX);
motion_feature_mask = zeros(h_motion,w_motion,'logical');

for i=1:N_motion_units
    motion_feature_mask(round(posX(i)), round(posY(i))) = true;
end

%% local mean subtraction and feature computations

gauss_size_pix = gauss_filter_size_static / pixel_resolution_static;
motion_center_size_pix = round(motion_feature_RF_center / pixel_resolution_motion);
motion_surround_size_pix = round(motion_feature_RF_surround/ pixel_resolution_motion);

% orientation
if features_static.horizontal_orientation || features_static.vertical_orientation
    OS_filt_h = fspecial('sobel');
    OS_filt_v = OS_filt_h';
end

if ~filter_gauss
    disk_filter = fspecial('disk', gauss_size_pix);
    disk_filter_motion_center = fspecial('disk', motion_center_size_pix);
    disk_filter_motion_surround = fspecial('disk', motion_surround_size_pix);
end

%initialize logical arrays for the features
feature_list_static = fieldnames(features_static);
for i=1:length(feature_list_static)
    feature_spikes_template.(feature_list_static{i}) = zeros(h_static,w_static,'logical');
end
feature_list_motion = fieldnames(features_motion);
for i=1:length(feature_list_motion)
    feature_spikes_template.(feature_list_motion{i}) = zeros(h_motion_feature,w_motion_feature,'logical');
end

feature_spikes = repmat(feature_spikes_template,frames);

%compute the features
tic;
for i=1:frames
    if filter_gauss
        f_filtered = imgaussfilt(frames_static(:,:,i),gauss_size_pix);
    else
        f_filtered = imfilter(frames_static(:,:,i), disk_filter);
    end
    f_background_subtracted = frames_static(:,:,i) - f_filtered;
    
    if features_static.local_contrast
        feature_spikes(i).local_contrast = abs(f_background_subtracted)>contrast_thres;
    end

    if features_static.horizontal_orientation
        feature_spikes(i).horizontal_orientation  = imfilter(f_background_subtracted,OS_filt_h,'replicate') > OS_thres;
    end

    if features_static.vertical_orientation
        feature_spikes(i).vertical_orientation  = imfilter(f_background_subtracted,OS_filt_v,'replicate') > OS_thres;
    end

    if features_motion.OMS
        if filter_gauss
            circ_filter_center = fspecial('gaussian',...
                [motion_center_size_pix*4, motion_center_size_pix*4],motion_center_size_pix);
            circ_filter_surround = fspecial('gaussian',...
                [motion_surround_size_pix*4, motion_surround_size_pix*4],motion_surround_size_pix);
        else
            circ_filter_center = fspecial('disk',motion_center_size_pix);
            circ_filter_surround = fspecial('disk',motion_surround_size_pix);
        end

        f_center_filtered = roifilt2(circ_filter_center, single(frames_diff(:,:,i)), motion_feature_mask);
        f_surround_filtered = roifilt2(circ_filter_surround, single(frames_diff(:,:,i)), motion_feature_mask);
        for j=1:N_motion_units
            x = round(posX(j));
            y = round(posY(j));
            feature_spikes(i).OMS(round(x/resize_factor_motion_feature), round(y/resize_factor_motion_feature)) = ...
                f_center_filtered(x,y) - ...
                motion_surround_weight * f_surround_filtered(x,y) > OMS_thres;
        end
    end
end
elapsed = toc;
fprintf('Feature extraction took %f seconds per frame\n', elapsed/frames);





