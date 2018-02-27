function [pIa, pIb] = match_sifts(Ia, Ib)

[fa,da] = vl_sift(im2single(Ia), 'PeakThresh', .01) ;
[fb,db] = vl_sift(im2single(Ib), 'PeakThresh', .01) ;

[matches, scores] = vl_ubcmatch(da,db) ;

[drop, perm] = sort(scores, 'descend') ;
matches = matches(:, perm) ;
scores  = scores(perm) ;

% matches = matches(:,1:100);
% scores = scores(:,1:100);

% do ransac
xa = fa(1,matches(1,:));
ya = fa(2,matches(1,:));

xb = fb(1,matches(2,:));
yb = fb(2,matches(2,:));


p1 = [ya; xa; ones(1,length(xa))];
p2 = [yb; xb; ones(1,length(xb))];

t = .005;  % Distance threshold for deciding outliers
%tic
% [F, inliers] = ransacfitfundmatrix(p1, p2, t);

[F, inliers] = ransacfithomography(p1, p2, t);
%toc

if 0
  if size(Ia,1) < size(Ib,1)
    figure(3);
    imshow( [[Ia ; zeros(size(Ib,1) - size(Ia,1), size(Ia,2))], Ib ]) ;
  else
    figure(3);
    imshow( [Ia , [Ib; zeros(size(Ia,1) - size(Ib,1), size(Ib,2))] ]) ;
  end
  
  xa = fa(1,matches(1,inliers)) ;
  xb = fb(1,matches(2,inliers)) + size(Ia,2) ;
  ya = fa(2,matches(1,inliers)) ;
  yb = fb(2,matches(2,inliers)) ;
  
  
  hold on ;
  h = line([xa ; xb], [ya ; yb]) ;
  set(h,'linewidth', 2, 'color', 'b') ;
end

pIa = [fa(1,matches(1,inliers)); fa(2,matches(1,inliers))];
pIb = [fb(1,matches(2,inliers)); fb(2,matches(2,inliers))];
