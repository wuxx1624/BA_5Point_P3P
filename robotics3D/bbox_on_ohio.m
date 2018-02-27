function [px, bbox] = bbox_on_ohio(coords)

% this function takes a lat-lon box, and returns the image patch on the
% ohio osip data

osip_org = [39.317233, -82.123561, 240];


% create an enu frame at the corner of the map
enu_R_ecef = ecef2lvRotationMatrix(osip_org(1,1) * pi/180, osip_org(1,2) * pi/180);

% the ecef coords of the corner
ecef_osip_org = lla2ecef(osip_org, 'WGS84')';

% the ps in enu
enu_p1 = enu_R_ecef*(coords(:,1) - ecef_osip_org);
enu_p2 = enu_R_ecef*(coords(:,2) - ecef_osip_org);
enu_p3 = enu_R_ecef*(coords(:,3) - ecef_osip_org);
enu_p4 = enu_R_ecef*(coords(:,4) - ecef_osip_org);


meters_to_feet = 3.280833333;
feet_to_pix = 1.0002;

px_p1 = enu_p1(1:2) * meters_to_feet * feet_to_pix;
px_p2 = enu_p2(1:2) * meters_to_feet * feet_to_pix;
px_p3 = enu_p3(1:2) * meters_to_feet * feet_to_pix;
px_p4 = enu_p4(1:2) * meters_to_feet * feet_to_pix;

% convert all pixel coordinates to the origin of the top left image in the
% array (i.e., s2075485)
px = [px_p1 px_p2 px_p3 px_p4];

px(2,:) = 10000 * ones(1,4) - px(2,:);

bbox.u = round([ min( px(1,:) ), max( px(1,:))]);
bbox.v = round([ min( px(2,:) ), max( px(2,:) )]);
% bbox.v = round([-max( [ px_p1(2) px_p2(2) px_p3(2) px_p4(2)] ) + 10000, -min( [ px_p1(2) px_p2(2) px_p3(2) px_p4(2)] ) + 10000]);


% make sure things are inbounds
for i = 1:2
  if bbox.u(i) < 1
    bbox.u(i) = 1;
  end
  if bbox.u(i) > 15000
    bbox.u(i) = 15000;
  end
  if bbox.v(i) < 1
    bbox.v(i) = 1;
  end
  if bbox.v(i) > 10000
    bbox.v(i) = 10000;
  end
end