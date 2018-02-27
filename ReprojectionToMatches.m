function [rmatch] = ReprojectionToMatches(P,R,t,CameraParams)


P1 = CameraParams.K*[R t];
Phom = [P(1:3,:); ones(1,size(P,2))];

rmatch = zeros(size(P,2), 2);

X1 = P1*Phom; rmatch(:, 1:2) = X1(1:2,:)'./repmat(X1(3,:)', [1 2]);



