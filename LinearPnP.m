function [R,t] = LinearPnP(P, w, CameraParams)

A = zeros(2*size(w,1), 12);

for i = 1:size(w,1)
    A(2*(i-1) + 1,1:4) = [P(1:3,i)' 1];
    A(2*(i-1) + 1,9:12) = -w(i,1)*[P(1:3,i)' 1];
    
    A(2*(i-1) + 2,5:8) = [P(1:3,i)' 1];
    A(2*(i-1) + 2,9:12) = -w(i,2)*[P(1:3,i)' 1];
end

[~, ~, V] = svd(A);

Pmat = [V(1:4,end)';V(5:8,end)';V(9:12,end)'];
[U,S,V] = svd(CameraParams.Kinv*Pmat(:,1:3));
R = U*V'; R = R/det(R);
t = CameraParams.Kinv*Pmat(:,4)/S(1,1);

