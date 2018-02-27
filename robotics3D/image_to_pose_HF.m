function [ecef_p_ned ned_p_body body_q_ned] = image_to_pose_HF(im_time)


load(strcat('c:/data/hf_data/truth/TrueTrajectory.mat'));


index = find(tgps < im_time, 1, 'last');

if (im_time - tgps(index)) / (tgps(index + 1) - tgps(index)) > 0.5
    index = index + 1;
end

    
ecef_p_ned = pe(1,:)';

ned_p_body = pn(index,:)';

ned_R_body = rotz(att(index,3)) * roty(att(index,2)) * rotx(att(index,1));

body_q_ned = rot2quat(ned_R_body');