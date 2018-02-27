function M = grobner_coefficient( EE )

syms x y z

% reshape the essential matrix to be 3*3
% E = zeros(3);
E = reshape((EE*[x;y;z;1]),[3,3]);

% determint constriant
M1(1,:) = (E(1,2)*E(2,3) - E(1,3)*E(2,2))* E(3,1)...
      + (E(1,3)*E(2,1) - E(1,1)*E(2,3))*E(3,2)...
      + (E(1,1)*E(2,2) - E(1,2)*E(2,1))*E(3,3);
  
A = 2*E*transpose(E)*E - trace(E*transpose(E))*E;
M1(2,:) = A(1,1);
M1(3,:) = A(1,2);
M1(4,:) = A(1,3);
M1(5,:) = A(2,1);
M1(6,:) = A(2,2);
M1(7,:) = A(2,3);
M1(8,:) = A(3,1);
M1(9,:) = A(3,2);
M1(10,:) = A(3,3);

M = zeros(10,20);
for i = 1:10
    M2 = coeffs(M1(i,:));
    M(i,:) = [M2(20) M2(19) M2(16) M2(10) M2(18) M2(15) M2(9) M2(13) M2(7) M2(4)...
    M2(17) M2(14) M2(8) M2(12) M2(6) M2(3) M2(11) M2(5) M2(2) M2(1)];
end