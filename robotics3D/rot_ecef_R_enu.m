function M = rot_ecef_R_enu(lat, lon)
% Construct the matrix that rotates Cartesian vectors from geocentric to
% local vertical.

% Copyright 2005 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2005/11/15 01:37:07 $

slat = sin(lat);
clat = cos(lat);
slon = sin(lon);
clon = cos(lon);

% M = [0 1 0; 1 0 0; 0 0 -1] * [    1      0        0   ; ...
%          0   slat  clat; ...
%          0  -clat  slat] ...
%   * [-slon   clon  0; ...
%      -clon  -slon  0; ...
%            0          0        1]; 
%        
%        
M = rotz(lon) * rotx(lat) * [0 0 1; 1 0 0; 0 1 0]';