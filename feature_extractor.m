function [] = feature_extractor(params)
% grab the parameters
%video - these are needed for saving the loaded video
pixel_resolution_static = params.video.pixel_resolution_static; %degrees
pixel_resolution_motion = params.video.pixel_resolution_motion; %degrees
video_fname = params.video.video_fname; %filename of the video (full path if not in current directory)
video_degrees_visual_angle = params.video.video_degrees_visual_angle;
contrast_max = params.video.contrast_max; % for 8 bit image, positive or negative contrast

video_params_hash = DataHash(params.video);

%feature extraction
gauss_filter_size_static = params.features.gauss_filter_size_static;
filter_gauss = params.features.filter_gauss;
OS_thres_frac = params.features.OS_thres_frac;
contrast_thres_frac = params.features.contrast_thres_frac;
features_static = params.features.features_static;
feature_spacing_motion = params.features.feature_spacing_motion;
motion_feature_RF_center = params.features.motion_feature_RF_center;
motion_feature_RF_surround = params.features.motion_feature_RF_surround;
motion_surround_weight = params.features.motion_surround_weight;
OMS_threshold_frac = params.features.OMS_threshold_frac;
loom_threshold_frac = params.features.loom_threshold_frac;
LD_contrast_thres_frac = params.features.LD_contrast_thres_frac;
DVS_contrast_thres_frac = params.features.DVS_contrast_thres_frac;
features_motion = params.features.features_motion;

all_params_hash = DataHash(params);

% %% parameters for static features
% pixel_resolution_static = .5; %degrees of visual angle
% gauss_filter_size_static = 1; %degrees of visual angle, s.d. of Gaussian or radius of disk
% filter_gauss = false; %true = gaussian, false = disk filter
% OS_thres_frac = .5; %fraction contrast
% contrast_thres_frac = .1;
%
% %which features
% features_static.local_contrast = true;
% features_static.horizontal_orientation = true;
% features_static.vertical_orientation = true;
%
% %% parameters for motion features
% pixel_resolution_motion = 0.2; %degrees
% feature_spacing_motion = 2; %degrees
% motion_feature_RF_center = 1; %degrees, could be fixed at spacing
% motion_feature_RF_surround = 5; %degrees
% motion_surround_weight = 2; %relative to center weight
% OMS_threshold_frac = .1; %fraction contrast
%
% %which features
% features_motion.OMS = true;
%
% %% Constants
% video_degrees_visual_angle = 120;
% contrast_max = 2^7; % for 8 bit image, positive or negative contrast.

%video_fname = 'GameRecording_short.mov';

%% Computed parameters
%contrast_max = 2^7;
OMS_thres = round(OMS_threshold_frac * contrast_max);
contrast_thres = round(contrast_thres_frac * contrast_max);
OS_thres = round(OS_thres_frac * contrast_max);
LD_contrast_thres = round(LD_contrast_thres_frac * contrast_max);
DVS_contrast_thres = round(DVS_contrast_thres_frac * contrast_max);

gauss_size_pix = gauss_filter_size_static / pixel_resolution_static;
motion_center_size_pix = round(motion_feature_RF_center / pixel_resolution_motion);
motion_surround_size_pix = round(motion_feature_RF_surround/ pixel_resolution_motion);

%% Load video
saved_frames_fname = sprintf('video_frames%s%s_vidFrames.mat', filesep, video_params_hash);
saved_video_params_fname = sprintf('video_frames%s%s_vidParams.mat', filesep, video_params_hash);

if exist(saved_frames_fname, 'file')
    load(saved_frames_fname)
    fprintf('Loaded saved movie franes from %s\n', saved_frames_fname);
else
    v = VideoReader(video_fname);
    degrees_per_pixel = video_degrees_visual_angle / v.Width;

    resize_factor_static = degrees_per_pixel / pixel_resolution_static;
    h_static = ceil(v.Height*resize_factor_static);
    w_static = ceil(v.Width*resize_factor_static);

    resize_factor_motion = degrees_per_pixel / pixel_resolution_motion;
    h_motion = ceil(v.Height*resize_factor_motion);
    w_motion = ceil(v.Width*resize_factor_motion);

    resize_factor_motion_feature = round(feature_spacing_motion / pixel_resolution_motion);
    h_motion_feature = ceil(h_motion/resize_factor_motion_feature);
    w_motion_feature = ceil(w_motion/resize_factor_motion_feature);

    frames_static = zeros(h_static,w_static,v.NumFrames,'uint8');
    frames_diff = zeros(h_motion,w_motion,v.NumFrames,'int8');

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

    fprintf('Saving movie frames to %s\n', saved_frames_fname);

    save(saved_frames_fname, 'frames_static', 'frames_diff', ...
        'resize_factor_static', 'resize_factor_motion', 'resize_factor_motion_feature', ...
        "h_motion_feature",'w_motion_feature','-v7.3');

    video_params = params.video;
    save(saved_video_params_fname, 'video_params');
end

[h_static, w_static, frames] = size(frames_static);
[h_motion, w_motion, ~] = size(frames_diff);

%% feature computations

saved_feature_spikes_fname = sprintf('feature_spikes%s%s_spikes.mat', filesep, all_params_hash);
saved_feature_params_fname = sprintf('feature_spikes%s%s_params.mat', filesep, all_params_hash);

if exist(saved_feature_spikes_fname, 'file')
    load(saved_feature_spikes_fname)
    fprintf('Loaded saved feature spikes from %s\n', saved_feature_spikes_fname);
else

    %% subunits for motion computations
    [posX, posY] = makeSubunitHexGrid([w_motion,h_motion], resize_factor_motion_feature);
    N_motion_units = length(posX);
    motion_feature_mask = zeros(h_motion,w_motion,'logical');

    for i=1:N_motion_units
        motion_feature_mask(round(posX(i)), round(posY(i))) = true;
    end

    %make binary mask for each motion unit center

    %gauss_size_pix
    motion_center_mask = false(h_motion, w_motion, N_motion_units);
    if features_motion.LD
        [X, Y] = meshgrid(1:w_motion,1:h_motion);
        for j=1:N_motion_units
            x = round(posX(j));
            y = round(posY(j));

            motion_center_mask(:,:,j) = sqrt((X - x).^2 + (Y - y).^2) < motion_center_size_pix;
        end
    end
    %% local mean subtraction and feature computations

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
        feature_spikes.(feature_list_static{i}) = zeros(h_static,w_static,frames,'logical');
    end
    feature_list_motion = fieldnames(features_motion);
    for i=1:length(feature_list_motion)
        feature_spikes.(feature_list_motion{i}) = zeros(h_motion_feature,w_motion_feature,frames,'logical');
    end

    %    feature_spike_table = table('Size',frames,)

    %feature_spikes = repmat(feature_spikes_template,frames);

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
            feature_spikes.local_contrast(:,:,i) = abs(f_background_subtracted)>contrast_thres;
        end

        if features_static.horizontal_orientation
            feature_spikes.horizontal_orientation(:,:,i)  = imfilter(f_background_subtracted,OS_filt_h,'replicate') > OS_thres;
        end

        if features_static.vertical_orientation
            feature_spikes.vertical_orientation(:,:,i)  = imfilter(f_background_subtracted,OS_filt_v,'replicate') > OS_thres;
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

            f_center_filtered = roifilt2(circ_filter_center, single(abs(frames_diff(:,:,i))), motion_feature_mask);
            f_surround_filtered = roifilt2(circ_filter_surround, single(abs(frames_diff(:,:,i))), motion_feature_mask);

            for j=1:N_motion_units
                x = round(posX(j));
                y = round(posY(j));
                feature_spikes.OMS(round(x/resize_factor_motion_feature), round(y/resize_factor_motion_feature), i) = ...
                    f_center_filtered(x,y) - ...
                    motion_surround_weight * f_surround_filtered(x,y) > OMS_thres;

            end
        end
        
        curFrameDiff = frames_diff(:,:,i);
        if features_motion.LD
            for j=1:N_motion_units
                x = round(posX(j));
                y = round(posY(j));
                positive_contrast_pix_count = sum(curFrameDiff(motion_center_mask(:,:,j)) > LD_contrast_thres,'all');
                negative_contrast_pix_count = sum(curFrameDiff(motion_center_mask(:,:,j)) < -LD_contrast_thres, 'all');
                %sum(motion_center_mask(:,:,j), 'all')
%                 if i==149
%                     keyboard;
%                 end

                feature_spikes.LD(round(x/resize_factor_motion_feature), round(y/resize_factor_motion_feature), i) = (negative_contrast_pix_count - positive_contrast_pix_count) / sum(motion_center_mask(:,:,j), 'all') > ...
                    loom_threshold_frac;

            end
        end
    end

    if features_motion.DVS
        feature_spikes.DVS = abs(frames_diff) > DVS_contrast_thres;
    end

    elapsed = toc;
    fprintf('Feature extraction took %f seconds per frame\n', elapsed/frames);
    
    fprintf('Saving feature spikes to %s\n', saved_feature_spikes_fname);
    save(saved_feature_spikes_fname,'feature_spikes', 'elapsed', '-v7.3')
    save(saved_feature_params_fname,'params');
end





