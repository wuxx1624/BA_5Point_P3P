function [px] = image_on_ohio(patch)

% this function takes a lat-lon box, and returns the image patch on the
% ohio osip data

osip_org = [39.317233, -82.123561, 240];


% create an enu frame at the corner of the map
enu_R_ecef = ecef2lvRotationMatrix(osip_org(1,1) * pi/180, osip_org(1,2) * pi/180);

% the ecef coords of the corner
ecef_osip_org = lla2ecef(osip_org, 'WGS84')';

% the image box in lat lon
p1 = [patch.lat(2), patch.lon(2), 240];
p2 = [patch.lat(2), patch.lon(1), 240];
p3 = [patch.lat(1), patch.lon(2), 240];
p4 = [patch.lat(1), patch.lon(1), 240];

% the point coordinates in ecef
ecef_p1 = lla2ecef(p1, 'WGS84')';
ecef_p2 = lla2ecef(p2, 'WGS84')';
ecef_p3 = lla2ecef(p3, 'WGS84')';
ecef_p4 = lla2ecef(p4, 'WGS84')';

% the ps in enu
enu_p1 = enu_R_ecef*(ecef_p1 - ecef_osip_org);
enu_p2 = enu_R_ecef*(ecef_p2 - ecef_osip_org);
enu_p3 = enu_R_ecef*(ecef_p3 - ecef_osip_org);
enu_p4 = enu_R_ecef*(ecef_p4 - ecef_osip_org);


meters_to_feet = 3.280833333;
feet_to_pix = 1.0002;

px_p1 = enu_p1(1:2) * meters_to_feet * feet_to_pix;
px_p2 = enu_p2(1:2) * meters_to_feet * feet_to_pix;
px_p3 = enu_p3(1:2) * meters_to_feet * feet_to_pix;
px_p4 = enu_p4(1:2) * meters_to_feet * feet_to_pix;

% convert all pixel coordinates to the origin of the top left image in the
% array (i.e., s2075485)
px.u = round([ min( [ px_p1(1) px_p2(1) px_p3(1) px_p4(1)] ), max( [ px_p1(1) px_p2(1) px_p3(1) px_p4(1)] )]);
px.v = round([-max( [ px_p1(2) px_p2(2) px_p3(2) px_p4(2)] ) + 10000, -min( [ px_p1(2) px_p2(2) px_p3(2) px_p4(2)] ) + 10000]);


% make sure things are inbounds
for i = 1:2
  if px.u(i) < 1
    px.u(i) = 1;
  end
  if px.u(i) > 15000
    px.u(i) = 15000;
  end
  if px.v(i) < 1
    px.v(i) = 1;
  end
  if px.v(i) > 10000
    px.v(i) = 10000;
  end
end