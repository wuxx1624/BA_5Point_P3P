function [pts2D, pts3D] = get_3D_image_coords(Ia)

[fa,da] = vl_sift(im2single(Ia)) ;


load Ground_sifts;

tic
[matches1, scores1] = vl_ubcmatch(da,dim1) ;
toc