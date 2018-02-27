% test search hf image region
clc; clear; close all;

vehicle = 'hf';
im_time = 406792.719820;

[ecef_p_ned ned_p_body body_q_ned] = image_to_pose_HF(im_time);
[patch] = image_on_earth( ecef_p_ned, ned_p_body, body_q_ned, vehicle );



clc; clear; close all;
% test search lf image region

im_time = 498775;
vehicle = 'lf';

[ecef_p_ned ned_p_body body_q_ned] = image_to_pose_LF(im_time);
[patch] = image_on_earth( ecef_p_ned, ned_p_body, body_q_ned, vehicle );