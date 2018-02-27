% This script using 2D matches to estimate position and orientation of new
% images

% TODO: another script to grasp the index automatically

clear all; close all;
run('/home/jahid/Fei/summer/vlfeat-0.9.20/toolbox/vl_setup.m');
addpath('./JacobianModule/');
addpath('./BundleAdjustment/');
addpath('./Drawings/');
addpath('./robotics3D/');

% Camera parameters:
% iphone 6

% fx = 1229; fy = 1153;
% cx = 360; cy = 640;

% K    = [fx 0 cx; ...
%         0 fy cy; ...
%         0 0 1];
% Kinv = [1/fx 0 -cx/fx; ...
%         0 1/fy -cy/fy; ...
%         0 0 1];

% vicon
fc = [254.997, 254.933];
cc = [326.69, 230.118];
kc = [0.925175, 0.0, 0.0, 0.0, 0.0];

fx = fc(1); fy = fc(2);
cx = cc(1); cy = cc(2);
K    = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];
% K = eye(3);
% Kinv = eye(3);

CameraParams = struct;
CameraParams.fx = fx; CameraParams.fy = fy; 
CameraParams.cx = cx; CameraParams.cy = cy; 
CameraParams.K = K; CameraParams.Kinv = Kinv;

% Control variables
% NumberOfImages = 5;
reloadImages = 1;
workingDirectory = './Set/';
reComputeMatches = 0;
pairWiseTriangulation = 0;
constructPoseGraphMatrix = 0;



matchesVisualization = 0;
reprojVisualization = 0;
multiPosesVisualization = 0;





%% Standard operations for obtaining pairwise matches and map builder
% Load images
fid = fopen('/home/jahid/documents/datasets/vicon_test/filter_result_1/image_names_1.txt','r');
bb = textscan(fid,'%s');
fclose(fid);

meaturePoses = [1;2;3];
NumberOfImages = length(meaturePoses);
I_rgb = cell(NumberOfImages, 1);
for i = 1:NumberOfImages
%     I_rgb{i} = imread([workingDirectory 'I' num2str(i) '.JPG']);
    I_rgb{i} = imread(['/home/jahid/documents',bb{1,1}{meaturePoses(i),1}]);
end

% Save or load feature extraction + desciptor
if reloadImages == 1
    featureExtracted = cell(NumberOfImages, 1);
    featureDescriptor = cell(NumberOfImages, 1);
    for i = 1:NumberOfImages
        if size(I_rgb{i},3)~=1
            I_rgb{i} = rgb2gray(I_rgb{i});
        end
        I_single = single(I_rgb{i});
        [featureExtracted{i},featureDescriptor{i}] = vl_sift(I_single, 'PeakThresh', 5, 'edgethresh', 5);
        figure;
        imshow(I_rgb{i});hold on
        plot(featureExtracted{i}(1,:),featureExtracted{i}(2,:), 'ro', 'markerSize', 5);
    end
    save([workingDirectory 'featureExtracted.mat'], 'featureExtracted');
    save([workingDirectory 'featureDescriptor.mat'], 'featureDescriptor');
else
    load([workingDirectory 'featureExtracted.mat']);
    load([workingDirectory 'featureDescriptor.mat']);
end


