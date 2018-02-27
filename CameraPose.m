function [R1, t1, R2, t2, R3, t3, R4, t4] = CameraPose(E)

[U,S,V]=svd(E);  
W=[0 -1 0;  
   1 0 0;  
   0 0 1  
   ];  

R1 = U*W*V';
R2 = U*W'*V';
if(det(R1)<0)
    R1 = -R1;
end
if(det(R2)<0)
    R2 = -R2;
end

t1 = U(:,3);
t2= - U(:,3);
t3 = t2;
t4 = t1;

R3 = R1;
R4 = R2;
