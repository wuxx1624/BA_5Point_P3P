function [C,p] = P3P_Solver(P_G_L,b_C)


if size(b_C,1) ~= 3
    b_C = [b_C;ones(1,size(b_C,2))];
end

%% Find the max ||(bi x bj)x(pi x pj)||
p12 = P_G_L(:,1)-P_G_L(:,2);
p23 = P_G_L(:,2)-P_G_L(:,3);
p31 = P_G_L(:,3)-P_G_L(:,1);

b_C(:,1) = b_C(:,1)/norm(b_C(:,1));
b_C(:,2) = b_C(:,2)/norm(b_C(:,2));
b_C(:,3) = b_C(:,3)/norm(b_C(:,3));
b1_b2 = cross(b_C(:,1),b_C(:,2));
b2_b3 = cross(b_C(:,2),b_C(:,3));
b3_b1 = cross(b_C(:,3),b_C(:,1));
%% Compute k1, k3'', ui, vipp
% k1,k2,k3,theta2
k1 = p12/norm(p12);
k3 = b1_b2/norm(b1_b2);
k2 = cross(k1,k3);
stheta2 = -k1'*k3;
ctheta2 = sqrt(1-stheta2^2);
k2 = k2/ctheta2;
% k3p
k3p = cross(k2,k1);
% cos(phi) sin(phi)
u1 = -p31;
u2 = p23;
% u1 = p31/norm(p31);
% u2 = p23/norm(p23);
csphi = [u1'*k2 u1'*k3p];
delta2 = csphi(1)^2 + csphi(2)^2;

% k3pp,vipp
% C_k2_theta2 = ctheta2*eye(3)-stheta2*skewsymm(k2)+(1-ctheta2)*(k2*k2');
v1 = -b3_b1;
v2 = b2_b3;
%% Compute f,g,A
% common terms used when checking the sign of sin(theta1')
k3b3 = k3'*b_C(:,3);
% f1i
u1k1 = (u1'*k1);
f14 = -u1k1*(k2'*b_C(:,1))*k3b3;
f15 = -u1k1*(k2'*v1);
f11 = (k2'*v1);
f12 = -(k2'*b_C(:,1))*k3b3;
f13 = (k3'*v1);
% f2i
u2k1 = (u2'*k1);
f24 = -u2k1*(k2'*b_C(:,2))*k3b3;
f25 = -u2k1*(k2'*v2);
f21 = (k2'*v2);
f22 = -(k2'*b_C(:,2))*k3b3;
f23 = (k3'*v2);
% gi
g1=f13*f22-f12*f23;
g2=f13*f25-f15*f23;
g3=f11*f23-f13*f21;
g4=f14*f23-f13*f24;
g5=f11*f22-f12*f21;
g6=f11*f25+f14*f22-f15*f21-f12*f24;
g7=f14*f25-f15*f24;
% Ai
A4=g5^2+g1^2+g3^2;
A3=2*(g5*g6+g1*g2+g3*g4);
A2=g6^2+2*g5*g7+g2^2+g4^2-(g1^2+g3^2)*delta2;
A1=2*g6*g7-2*(g1*g2+g3*g4)*delta2;
A0=g7^2-(g2^2+g4^2)*delta2;

if isnan(A0) || A0 == Inf || isnan(A1) || A1 == Inf || isnan(A2) || A2 == Inf || isnan(A3) || A3 == Inf || isnan(A4) || A4 == Inf
    C=[];
    p=[];
else
%% Solve cos(theta1')
ss = roots([A4 A3 A2 A1 A0]);

if isempty(ss)
    C=[];
    p=[];
else
%% Confirm cos(theta1')
% % Ground truth
% u=k1;
% v=k3;
% [s1,~]=find_theta(C_G_C,u,v);
% % show
% ss
% cos(s1(1)-atan2(csphi(2),csphi(1)))
%% Get C
% compute C for all real cos(theta1')
C=[];
p = [];
k=1;
Binvb = b_C(:,3)/k3b3;
real_num = 0;


for i=1:length(ss)
   % check if it is a real solution
   if isreal(ss(i))
       ctheta1p(k) = ss(i);
       real_num = real_num + 1;
       k = k+1;
   else
       continue;
   end
end

% if real_num == 0
%     for i = 1:4
%         [~,idx_ss] = min(imag(ss(i)));
%     end
%     for i = 1 : length(idx_ss)
%         ctheta1p(i) = ss(i);
%     end
% end

if real_num == 0
    [~,idx_ss] = min(imag(ss));
    for i = 1 : length(idx_ss)
        ctheta1p(i) = ss(i);
    end
end

for k=1:length(ctheta1p)
   if k3b3>0
       stheta1p = sqrt(delta2-ctheta1p(k)^2);
   else
       stheta1p = -sqrt(delta2-ctheta1p(k)^2);
   end
   % cos(theta3) sin(theta3)
   cstheta3 = [g1*ctheta1p(k)+g2;g3*ctheta1p(k)+g4]*(stheta1p/(g5*ctheta1p(k)^2+g6*ctheta1p(k)+g7));
   % cos(theta1) sin(theta1)
   ctheta1 = (ctheta1p(k)*csphi(1)-stheta1p*csphi(2))/delta2;
   stheta1 = (stheta1p*csphi(1)+ctheta1p(k)*csphi(2))/delta2;
   % compute C
   ctheta3=cstheta3(1);
   stheta3=cstheta3(2);
   mc = [ctheta3 0 -stheta3;
       stheta1*stheta3 ctheta1 stheta1*ctheta3;
       ctheta1*stheta3 -stheta1 ctheta1*ctheta3];
   fc = [k1 k3p k2];
   bc = [cross(k3,k2)';k3';k2'];
   C(:,:,k)=fc*mc*bc;
   % Get p
   p(:,k) = [P_G_L(:,3)-C(:,:,k)*stheta1p*Binvb];
   
   % change reference order
   C(:,:,k) = C(:,:,k)';
   p(:,k) = -C(:,:,k)*p(:,k);
end
end
end
