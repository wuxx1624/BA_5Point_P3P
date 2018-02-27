function [ecef_p_ned ned_p_body body_q_ned] = image_to_pose_LF(im_time)


 load(strcat('../../lf_data/truth/TrueTrajectory.mat'));
 % load /home/joel/data/lf_data/truth/TrueTrajectory.mat;


% this file gives
% ecef_p_ned, the origin of the ned frame in ecef
% ned_p_body, the trajectory of the body in ned
% body_q_ned, the orientation of the ned in the body frame
% times, the gps time of the estimates


index = find(tgps < im_time, 1, 'last');


if (im_time - tgps(index)) / (tgps(index + 1) - tgps(index)) > 0.5
    index = index + 1;
end

ecef_p_ned = ecef_p_ned';

ned_p_body = ned_p_body(:,index); 

body_q_ned = body_q_ned(:,index);