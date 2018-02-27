% XYZ2GRAV - computes the gravitational acceleration vector at a specified
%            ECF location using the JGM2 gravitational ellipsoid only.
%            Higher-order gravity terms (the "gravity anomaly") are
%            ingnored. Only the pure ellipsoid is used.
% 
% [gx,gy,gz] = XYZ2GRAV(x,y,z)
% 
%  x, y, z = Earth-Centered-Fixed (ECF) cartesian coordinates (meters)
% gx,gy,gz = components of gravitational acceleration in ECF (meters/sec^2)
% 
% NOTES:  (1) x,y,z may be scalars, vectors, or matrices but must have
%             the same size and shape
%         (2) gx,gy,gz will have the same size/shape as x,y,z
%         (3) Only the purely ellipsoidal portion of the JGM-2 gravity
%             model is used here
%         (4) Points inside the earth will return gravity vectors, but
%             will be higher in magnitude than physical reality since
%             mass "above" the point is not removed from consideration;
%             i.e., the ellipsoidal gravitational field surrounds a point at
%             the center of the earth. The strength of the field increases
%             as one gets closer to that central point.
%         (5) If all inputs are zero, then NaNs are returned.
%         (6) Source of initial formulas:
%             www.colorado.edu/ASEN/asen3200/labs/ASEN3200_LabO3_2005.pdf
%         (7) The gx & gy formulas were written in a more explicit form
%             to prevent divergence when the x-coordinate is zero.
%         (8) The calculation was accelerated by precomputing some values.
%         (9) No warranty; use at your own risk.
%        (10) Version 1.0 - Initial writing. Michael Kleder, August 2005
%
% % EXAMPLE:
% R=6378137;
% [x,y,z]=sphere(50);
% [x,y,z]=deal(x*R,y*R,z*R);
% [gx,gy,gz]=xyz2grav(x,y,z);
% gm = sqrt(gx.^2+gy.^2+gz.^2);
% figure
% surf(x,y,z,gm)
% axis equal
% colorbar
% title('Gravitational acceleration on a uniform sphere (m/s^2)')
% disp(['When standing at a pole, one experiences slightly MORE gravitational ' ...
%     'acceleration, but that is because one is closer to the center of ' ...
%     'the earth (smaller radius); however, the provided plot is of ' ...
%     'gravity on a uniform sphere -- equal distances -- so that gravity ' ...
%     'is greater along the equatorial plane.'])

function [gx,gy,gz] = xyz2grav(x,y,z)
J2 = .00108263;
mu = 3.986004418e14; % m^3/s^2
R = 6378137; % m
r = sqrt(x.^2+y.^2+z.^2);
sub1 = 1.5*J2.*(R./r).^2;
sub2 = 5.*z.^2./r.^2;
sub3 = -mu./r.^3;
sub4 = sub3.*(1-sub1.*(sub2-1));
gx = x .* sub4;
gy = y .* sub4;
gz = z .* sub3.*(1-sub1.*(sub2-3));
return
