function [height] = dted_height(lla)

lat_deg = floor(abs(lla(1)));
lat_dec = rem(abs(lla(1)), 1);

lon_deg = floor(abs(lla(2)));
lon_dec = rem(abs(lla(2)), 1);

dir1 = strcat('W0', num2str(lon_deg));

file1 = strcat('N', num2str(lat_deg), '.DT1');

dted_file = strcat('c:/data/DTED1/', dir1, '/' , file1);

[Z, refvec, UHL, DSI, ACC] = dted(dted_file);


height = Z( round(size(Z,1) * lat_dec),  round(size(Z,2) * lon_dec)) + geoidheight(lla(1), lla(2), 'None');