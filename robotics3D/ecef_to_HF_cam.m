function [px, bbox] = ecef_to_HF_cam(ecef_pt, ecef_p_hf_ned, ned_p_hf_body, body_q_hf_ned, calib_switch)

px = zeros(2,4);

% compute the corresponding lat long alt
lla_p_hf_ned = ecef2lla(ecef_p_hf_ned');

% rotation from ecef to the east north up frame located at pe(1)
enu_R_ecef = ecef2lvRotationMatrix(lla_p_hf_ned(1,1) * pi/180, lla_p_hf_ned(1,2) * pi/180);

% rotation from east north up to north east down, this is useful
ned_R_enu = [0 1 0; 1 0 0; 0 0 -1];

ned_R_ecef = ned_R_enu * enu_R_ecef;

hf_R_ned = quat2rot(body_q_hf_ned);

%hf_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752;
%            0.004644596955214  0.986541150073639    0.16344719921;
%            0.000356118382038 -0.163450583641475    0.98655145830]

if strcmp(calib_switch, 'their_calib')
    hf_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752128807;
                0.00464459695521431  0.986541150073639    0.163447199212791;
                0.000356118382038752 -0.163450583641475    0.98655145830];
    hf_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752128807;
                0.00464459695521431 0.986541150073639 0.163447199212791;
                0.000356118382038752 -0.163450583641475 0.986551458306641] * rotz(pi/2);
   % [U S V] = svd(hf_R_cam);
    
%    hf_R_cam = U*V* rotz(pi/2);
            

     %       ;
    
    
    fc =  [14064.38769;
           14064.38769];
       
    cc = [2436; 1624];% + [-465 ; -1281];
elseif strcmp(calib_switch, 'our_calib')
    % my calibration
    
    load '../Joel-hf-cam-calib-extrinsic/body_R_cam.mat';
    hf_R_cam = quat2rot(q_est);
    
    rpy = rot2rpy(hf_R_cam);
    
%    droll = 0; -1.95;
%    dpitch = 0; 2.5;
%    dyaw = 0; -5;
    
    
    % hf_R_cam = rotz(rpy(3) + dyaw * pi/180) * roty(rpy(2) + dpitch*pi/180) * rotx(rpy(1) + droll*pi/180);
    
    
    hf_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752;
        0.004644596955214  0.986541150073639    0.16344719921;
        0.000356118382038 -0.163450583641475    0.98655145830]* rotz(pi/2);
    
    dyaw = -.05;
    dpitch = 0.87;
    droll = -3;
    
    hf_R_cam = hf_R_cam * rotz(dyaw * pi/2) * roty( dpitch * pi/180 ) * rotx(droll * pi/180);
    
    load '../Joel-hf-cam-calib/Calib.mat' 
    

%     
%         fc =  1.1 * [14064.38769;
%         14064.38769];
%     
%     cc = [2436; 1624];
end

hf_p_cam = zeros(3,1);

for i=1:4
    ned_pt = ned_R_ecef * (ecef_pt(:,i) - ecef_p_hf_ned);
    
    hf_pt = hf_R_ned * ( ned_pt - ned_p_hf_body );
    
    cam_pt = hf_R_cam' * ( hf_pt - hf_p_cam);
    
    px(:,i) = [ fc(1) * cam_pt(1) / cam_pt(3) + cc(1);
        fc(2) * cam_pt(2) / cam_pt(3) + cc(2)];
    
end

bbox.u = [min(px(1,:)) max(px(1,:))];
bbox.v = [min(px(2,:)) max(px(2,:))];