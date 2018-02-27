function [I2pts3D] = get_3Dpts_ohio(px, I2pts)

foot_to_meter = 1 / 3.280833333;

load './data/X3D.mat';

X3D = X3D(px.v(1):px.v(2), px.u(1):px.u(2));

I2pts3D = zeros(3, length(I2pts));

for i = 1:length(I2pts)
  I2pts3D(1,i) = X3D( round(I2pts(2,i)), round(I2pts(1,i)));
end

clear X3D;

% for y we save the transpose distances matrix since it's more efficient
load './data/Y3Dt.mat';
% Y3D = Y3D(px.v(1):px.v(2), px.u(1):px.u(2));
Y3D = Y3Dt(px.u(1):px.u(2), px.v(1):px.v(2))'; 
for i = 1:length(I2pts)
  I2pts3D(2,i) = Y3D( round(I2pts(2,i)), round(I2pts(1,i)));
end
clear Y3Dt;

load './data/Z3D.mat';
Z3D = Z3D(px.v(1):px.v(2), px.u(1):px.u(2));
for i = 1:length(I2pts)
  I2pts3D(3,i) = Z3D( round(I2pts(2,i)), round(I2pts(1,i)));
end
clear Z3D;


I2pts3D = I2pts3D * foot_to_meter;

% center the data, this will help with conditioning in the camera
% calibration
I2pts3D = I2pts3D - repmat(mean(I2pts3D,2), 1, length(I2pts3D));