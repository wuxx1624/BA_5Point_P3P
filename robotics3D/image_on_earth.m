function [patch, lla_p_body, ecef_pt] = image_on_earth( ecef_p_ned, ned_p_body, body_q_ned, vehicle , calib)

addpath(genpath( '../../TOOLBOX_calib/') );

% inputs:
% pe is the position of the nav frame origin in ecef
% pn is the vehicle position in the nav frame (i.e., north-east-down frame), specifically nav_p_body
% q is the vehicle orientation in the nav frame, specifically body_q_nav
%
% outputs:
% patch.lat = [min lat, max lat]
% patch.lon = [min lon, max lon]

% This code works by
% 1. creating a east-north-up frame directly below the vehicle on the
% earth's surface
% 2. taking four unit vectors, corresponding to the four corners of the
% image
% 3. ray-tracing the unit vectors from the camera to the intersection with
% the earth, assuming that it is locally flat
% 4. computing the bounding box in the local east-north-up frame
% 5. converting those coordinates back to lat-lon-alt

if strcmp(vehicle, 'hf')
  % fov is H x V = 20 deg x 13 deg
  % alpha = 19.65/2 * pi/180;
  % beta = 13.17/2 * pi/180;
  
  alpha = 13.17/2 * pi/180; % based on GFE calibration, vertical fov
  beta = 19.65/2 * pi/180; % based on GFE calibration, horizontal fov
  
  
  if strcmp(calib, 'their_calib')
  body_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752;
                0.004644596955214  0.986541150073639    0.16344719921;
                0.000356118382038 -0.163450583641475    0.98655145830] * rotz(pi/2);
            
  elseif strcmp(calib, 'our_calib')
    load '../Joel-hf-cam-calib-extrinsic/body_R_cam.mat';
    body_R_cam = quat2rot(q_est);
    
      body_R_cam = [0.999989150390554 -0.00452392734728083 -0.00111048752;
                0.004644596955214  0.986541150073639    0.16344719921;
                0.000356118382038 -0.163450583641475    0.98655145830] * rotz(pi/2);
  end
  
  fc = [14064.38769; 14064.38769];
  cc = [2436; 1624];
  kc = zeros(5,1);
  alpha_c = 0;
  
      pt1 = [1;1];
    pt2 = [4872; 1];
    pt3 = [4872 ; 3248];
    pt4 = [1 ; 3248];
    
    pts = [pt1 pt2 pt3 pt4];
            
% my calibration
% body_R_cam = rotz(86.007091964807756 * pi/180) * roty( 23.709790784135908 * pi/180) * rotx(-1.003118193730841 * pi/180);
 

body_p_cam = zeros(3,1);

pt_n = normalize(pts, fc, cc, kc, alpha_c);
  
  pt_n = [pt_n ; ones(1,4)];
  
  pts = pt_n ./ repmat(sqrt(sum(pt_n.*pt_n,1)), 3,1);
  
elseif strcmp(vehicle, 'lf')
  
    pt1 = [1;1];
    pt2 = [1360; 1];
    pt3 = [1360 ; 1024];
    pt4 = [1 ; 1024];
    
    pts = [pt1 pt2 pt3 pt4];
    
    if strcmp(calib, 'their_calib')
    % fov is H x V = 56.18 deg x 43.79 deg
%   alpha = 43.79/2 * pi/180;  % based on the gfe calibration
%   beta = 56.18/2 * pi/180; % based on the gfe calibration
  
  pp = [741.583; 528.537];
  foc = 1274.034;
  
  pts = [(pts - repmat(pp, 1,4))/foc ; ones(1,4)];
  
  pts = pts ./ repmat(sqrt(sum(pts.*pts,1)),3,1);
  
  body_R_cam = rotz(0.224467 * pi/180) * roty( -4.399699 * pi/180) * rotx( 1.670560 * pi/180) * rotz(pi/2);
  
elseif strcmp(calib, 'our_calib')
  % fov is H x V = 46.70 deg x 36.01 deg
  alpha = 36.01/2 * pi/180;  % based on my calibration, vertical fov
  beta = 46.70/2 * pi/180; % based on my calibration, horizontal fov
  
  
  cc = [772.4114; 484.11];
  fc = .85 * [1584.15; 1581.47];
  kc = [-.18972 0 0 0 0 0];
  
  alpha_c = 0;
  
  pt_n = normalize(pts, fc, cc, kc, alpha_c);
  
  pt_n = [pt_n ; ones(1,4)];
  
  pts = pt_n ./ repmat(sqrt(sum(pt_n.*pt_n,1)), 3,1);
  
  % this is my calibration
%  body_R_cam = rotz(89.694380810884795 * pi/180) * roty( 1.117100762053628 * pi/180) * rotx( -3.833994025516499 * pi/180);
%  body_R_cam = rotz(89.86 * pi/180) * roty( 1.42 * pi/180) * rotx( -3.48 * pi/180);
  
  load '../Joel-lf-cam-calib-extrinsic/body_R_cam.mat';
  body_R_cam = quat2rot(q_est);
  
  body_R_cam = rotz(89.7213 * pi/180) * roty( 1.1499 * pi/180) * rotx( -3.7401 * pi/180);
  
  end
  
 % body_p_cam = [1.212; -0.368; 2.509] - [6.087; 0.448; 17.62];
  
end
% z_axis = [0 0 1]';
% corner1 = roty(beta) * rotx(alpha) * z_axis;
% corner2 = roty(beta) * rotx(-alpha) * z_axis;
% corner3 = roty(-beta) * rotx(alpha) * z_axis;
% corner4 = roty(-beta) * rotx(-alpha) * z_axis;
% 
% 
% points = [corner1 corner2 corner3 corner4];

if strcmp(calib, 'their_calib')
    points = pts;
elseif strcmp(calib, 'our_calib')
    points = pts;
end

ground_height = 200; % this is approximate

% compute the lla of the nav frame origin
lla = ecef2lla(ecef_p_ned');

% rotation from ecef to the east north up frame located at pe
enu_R_ecef = ecef2lvRotationMatrix(lla(1,1) * pi/180, lla(1,2) * pi/180);

ned_R_enu = [0 1 0; 1 0 0; 0 0 -1];

ned_R_ecef = ned_R_enu * enu_R_ecef;

% compute the vehicle position in ecef
ecef_p_body = ecef_p_ned + ned_R_ecef' * ned_p_body;

% compute the vehicle position in lla
lla_p_body = ecef2lla(ecef_p_body');

% create an enu frame directly under the vehicle on the surface
enu2_R_ecef = ecef2lvRotationMatrix(lla_p_body(1,1) * pi/180, lla_p_body(1,2) * pi/180);

lla_pt = zeros(4,3);
ecef_pt = zeros(3,4);

for i = 1:4
  body_unit_vector = body_R_cam * points(:,i); % note that the camera and the body frame are coincident
  
  enu_unit_vector = enu2_R_ecef * ned_R_ecef' * quat2rot(body_q_ned)' * body_unit_vector;
  
  theta = acos(-enu_unit_vector(3));
  
  d = tan(theta) * (lla_p_body(1,3)-ground_height); % the ground is approx 200 m above the ellipsoid around our nav region
  
  enu_xy_unit_vector = enu_unit_vector(1:2) / norm(enu_unit_vector(1:2));
    
  enu2_pt = [ d * enu_xy_unit_vector ; 
              0];
  ecef_pt(:,i) = ecef_p_body + enu2_R_ecef' * (enu2_pt - [0 0 (lla_p_body(1,3)-ground_height)]');
  
  lla_pt(i,:) = ecef2lla(ecef_pt(:,i)'); % find the point's coordinates in lat-lon-alt
  lla_pt(i,3) = 0; % project the point back to the earth's surface
                   
end

patch.lat = [min(lla_pt(:,1)) max(lla_pt(:,1))];
patch.lon = [min(lla_pt(:,2)) max(lla_pt(:,2))];