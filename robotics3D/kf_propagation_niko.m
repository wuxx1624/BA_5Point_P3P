function [xek_1k,Pk_1k,Fk,Qdk]=kf_propagation_niko(xekk,Pkk,GYRO_previous,GYRO_current,ACCEL_previous,ACCEL_current,DT,g,sigma_r,sigma_w,sigma_g,sigma_v)

% INPUTS: 
% 1. previous step values (or initial) of STATE ESTIMATES
% 2. previous step values (or initial) of COVARIANCE ESTIMATES 
% 3. IMU sensor measurements( GYRO & ACCELS) - PREVIOUS & CURRENT TIME STEP
% These are expressed in the local frame
% ACCEL vector contains GRAVITATIONAL ACCELERATION vector & BIASES
% 4. noise sigma|s for IMU sensors -> white + bias 
% 5. DT: time step
% 6. gravitational acceleration (in global coordinates)
% 7. functions needed for integration processes

% OUTPUTS:
% 1. next step values of STATE ESTIMATES
% 2. next step values of COVARIANCE ESTIMATES

%-------------------------------------------------------------------
% xekk - STATE ESTIMATES (16 state estimates  propagated - 15 independent)
%
% from GLOBAL to LOCAL - ORIENTATION
% 1 q1
% 2 q2
% 3 q3
% 4 q4
%
% SENSOR BIAS - GYRO
% 5 gyro_bias a
% 6 gyro_bias b
% 7 gyro_bias c
%
% wrt GLOBAL - VELOCITY
% 8 v1 
% 9 v2
% 10 v3
%
% SENSOR BIAS - ACCEL
% 11 accel_bias a
% 12 accel_bias b
% 13 accel_bias c
%
% wrt GLOBAL - POSITION
% 14 r1
% 15 r2
% 16 r3
%
%-----------------------------------------------------------------------
%
% Pkk - COVARIANCE ESTIMATES
%
% FROM GLOBAL TO LOCAL - ORIENTATION ERROR
% 1 dq1 *2 = d(theta1)
% 2 dq2 *2 = d(theta2)
% 3 dq3 *2 = d(theta3)
%
% SENSOR BIAS - GYRO ERROR
% 4 D(gyro_bias a)
% 5 D(gyro_bias b)
% 6 D(gyro_bias c)
%
% wrt GLOBAL VELOCITY ERROR
% 7 D(v1)
% 8 D(v2)
% 9 D(v3)
%
% SENSOR BIAS - ACCEL ERROR
% 10 D(accel_bias a)
% 11 D(accel_bias b)
% 12 D(accel_bias c)
%
% wrt GLOBAL - POSITION ERROR
% 13 D(r1)
% 14 D(r2) 
% 15 D(r3)
%
%-----------------------------------------------------------------------
%///////////////////////////////////////////////////////////////////////
%-----------------------------------------------------------------------
%
%    S T A T E    E S T I M A T E S    P R O P A G A T I O N 
%
%-----------------------------------------------------------------------
%///////////////////////////////////////////////////////////////////////
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% GYRO BIAS PROPAGATION (PREDICTION)
%-----------------------------------------------------------------------
% Previous bias estimate - tk
gyro_bias_kk = xekk(5:7,1);

% Current bias estimate - tk+1
gyro_bias_k_1k = gyro_bias_kk;

% propagated state
xek_1k(5:7,1) = gyro_bias_k_1k;

%-----------------------------------------------------------------------
% ORIENTATION PROPAGATION (PREDICTION)
%-----------------------------------------------------------------------
% Previous UNBIASED rotational velocity estimate
w_kk = GYRO_previous - gyro_bias_kk;

% Current UNBIASED rotational velocity estimate
w_k_1k = GYRO_current - gyro_bias_k_1k;

% average rotational velocity estimate
w_avg = 0.5*(w_kk + w_k_1k);

% /// CHANGE begins here (Nikolas, July 9, 2005) (moved this portion up
% from further down to also use it during vel/pos integration)
% norm of avg. rotational velocity
norm_w_avg = sqrt(w_avg' * w_avg);
if norm_w_avg < 0.01
 %   disp('SMALL NORM OF OMEGA!')
end

% negative cross product matrix derived using GYRO measurements
[NCP_w_avg]=neg_cross_product(w_avg);

% /// CHANGE ends here


% Skew symmetric matrices calculations
[Omega_a]=Omega_calculation(w_avg);
[Omega_kk]=Omega_calculation(w_kk);
[Omega_k_1k]=Omega_calculation(w_k_1k);

% Previous orientation (QUATERNIONS) estimate
qkk= xekk(1:4,1);

% orientation prediction
% Note: Equations from page-123-Thesis
qk_1k = [ expm(Omega_a * DT/2) + ( Omega_k_1k * Omega_kk - Omega_kk * Omega_k_1k ) * ( DT^2)/48 ] * qkk;

% // CHANGE - addition of the following lines to keep the 4th quaternion element always positive//
if ~(abs(qk_1k(4,1))==0)
    qk_1k=sign(qk_1k(4,1))*qk_1k;
end

% quaternion normalization
qk_1k = qk_1k/(sqrt(qk_1k'*qk_1k));

% propagated state
xek_1k(1:4,1) = qk_1k;

%-----------------------------------------------------------------------
% ACCEL BIAS PROPAGATION (PREDICTION)
%-----------------------------------------------------------------------
% Previous bias estimate - tk
accel_bias_kk = xekk(11:13,1);

% Current bias estimate - tk+1
accel_bias_k_1k = accel_bias_kk;

% propagated state
xek_1k(11:13,1) = accel_bias_k_1k;

%-----------------------------------------------------------------------
% VELOCITY PROPAGATION (PREDICTION)
% LOCAL ACCELERATION INTEGRATION -> GLOBAL VELOCITY
%-----------------------------------------------------------------------
% Previous velocity estimate - tk
vkk = xekk(8:10,1);

% Previous UNBIASED acceleration estimate - tk
dot_v_kk = ACCEL_previous - accel_bias_kk;

% Current UNBIASED acceleration estimate - tk+1
dot_v_k_1k = ACCEL_current - accel_bias_k_1k;

% /// CHANGE begins here (Nikolas, July 9, 2005) (added this quantity)
% Average UNBIASED acceleration estimate
dot_v_avg = 0.5 * (dot_v_kk + dot_v_k_1k);

% /// CHANGE begins here (Nikolas, July 9, 2005) (this portion not needed)
% % inverse quaternion at tk
% [inv_qkk]=quaternion_inv(qkk);
% 
% % inverse quaternion at tk+1
% [inv_qk_1k]=quaternion_inv(qk_1k);
% 
% % quaternion-like Global gravitational acceleration
% G_g_q=[g;
%        0];
% 
% % Previous gravitational acceleration (local) - tk
% [L_gkk_q]=quaternion_mul(qkk,G_g_q);
% [L_gkk_q]=quaternion_mul(L_gkk_q,inv_qkk);
% 
% % Current gravitational acceleration (local) - tk+1
% [L_gk_1k_q]=quaternion_mul(qk_1k,G_g_q);
% [L_gk_1k_q]=quaternion_mul(L_gk_1k_q,inv_qk_1k);
% 
% % Previous UNBIASED & GRAVITATION-FREE acceleration estimate - tk
% dot_v_kk = dot_v_kk - L_gkk_q(1:3,1);
% 
% % Current UNBIASED & GRAVITATION-FREE acceleration estimate - tk
% dot_v_k_1k = dot_v_k_1k - L_gk_1k_q(1:3,1);
% 
% % Average UNBIASED & GRAVITATION-FREE acceleration estimate
% dot_v_avg = 0.5 * (dot_v_kk + dot_v_k_1k);
% % quaternion-like average UNBIASED & GRAVITATION-FREE acceleration estimate
% dot_v_avg_q=[dot_v_avg;
% 	     0];

% /// CHANGE ends here


% Skew symmetric matrix calculation: w_avg is the previously
% calculated average rotational velocity during [tk,tk+1]
[para_Omega_a]=para_Omega_calculation(w_avg);

% INTEGRATION OF (LOCAL) ACCELERATION -> (GLOBAL) VELOCITY
%----------------------------
% /// CHANGE begins here (Nikolas, July 9, 2005) (New Integrator)
% % FIRST TERM
% % Dv1= q^(-1) (x) a (x) q * DT
% [Dv1]=quaternion_mul(inv_qkk,dot_v_avg_q);
% [Dv1]=quaternion_mul(Dv1,qkk);
% Dv1=Dv1 * DT;
% 
% % SECOND TERM
% % q^(-1) (x) a (x) (Omega(w)*q) * (DT^2)/4
% [Dv2]=quaternion_mul(inv_qkk,dot_v_avg_q);
% [Dv2]=quaternion_mul(Dv2,Omega_a*qkk);
% Dv2 = 0.25 * (DT^2) * Dv2;
% 
% % THIRD TERM
% % (para_Omega(w)*q^(-1)) (x) a (x) q * (DT^2)/4
% [Dv3]=quaternion_mul(para_Omega_a*inv_qkk,dot_v_avg_q);
% [Dv3]=quaternion_mul(Dv3,qkk);
% Dv3 = 0.25 * (DT^2) * Dv3;
% 
% % TOTAL VELOCITY
% vk_1k = vkk + Dv1(1:3,1) + Dv2(1:3,1) + Dv3(1:3,1);

% For this integrator, cf. Quaternion Techreport, section 4.3 (Integration
% of Position and Velocity)
% FIRST TERM
if norm_w_avg < 0.1
    %small omega solution
    Dv1 = quaternions_to_rot_matrix(qk_1k)'...
        * [ DT*eye(3)  +  NCP_w_avg * 1/2*DT^2  +  NCP_w_avg * NCP_w_avg  * 1/6*DT^3 ]...
        * dot_v_avg;
else 
    %normal solution
    Dv1 = quaternions_to_rot_matrix(qk_1k)'...
        * [ DT*eye(3)  +  NCP_w_avg * 1/norm_w_avg^2 * ( 1-cos(norm_w_avg*DT) ) ...
        +  NCP_w_avg * NCP_w_avg  * 1/norm_w_avg^3 * ( norm_w_avg*DT - sin(norm_w_avg*DT) ) ]...
        * dot_v_avg;
end
    
% SECOND TERM (acceleration due to gravity)
Dv2 = +DT * g;

% TOTAL VELOCITY
vk_1k = vkk + Dv1 + Dv2;

% /// CHANGE ends here


% estimated state propagation tk+1
xek_1k(8:10,1) = vk_1k;

%-----------------------------------------------------------------------
% POSITION PROPAGATION
%-----------------------------------------------------------------------
% --- Useful quantities---
% average estimated velocity - during [tk,tk+1]
v_avg=0.5*(xekk(8:10,1) + xek_1k(8:10,1));
% quaternion-like average estimated velocity - during [tk,tk+1]
v_avg_q=[v_avg;
	   0];

% previous position estimate (global)
rkk = xekk(14:16,1);

% INTEGRATION OF GLOBAL VELOCITY & LOCAL ACCELERATION -> GLOBAL POSITION
%-----------------------------------
% /// CHANGE begins here (Nikolas, July 9, 2005) (New Integrator)
% % FIRST TERM
% % Dr1 = v * DT
% Dr1 = v_avg_q * DT;
% 
% % SECOND TERM
% % Dr2 = q^(-1) (x) a (x) q * (DT^2)/2
% [Dr2]=quaternion_mul(inv_qkk,dot_v_avg_q);
% [Dr2]=quaternion_mul(Dr2,qkk);
% Dr2 = 0.5 * (DT^2) * Dr2;
% 
% % THIRD TERM - A
% % Dr3a  = (para_Omega(w)*q^(-1)) (x) q (x) v * (DT^2)/4
% [Dr3a]=quaternion_mul(para_Omega_a*inv_qkk,qkk);
% [Dr3a]=quaternion_mul(Dr3a,v_avg_q);
% Dr3a = 0.25 * (DT^2) * Dr3a;
% 
% % THIRD TERM - B
% % Dr3b = v (x) q (x) (Omega(w)*q) * (DT^2)/4
% [Dr3b]=quaternion_mul(v_avg_q,qkk);
% [Dr3b]=quaternion_mul(Dr3b,Omega_a*qkk);
% Dr3b = 0.25 * (DT^2) * Dr3b;
% 
% rk_1k=rkk + Dr1(1:3,1) + Dr2(1:3,1);
% % *** Keep these terms out for now + Dr3a(1:3,1) + Dr3b(1:3,1);

% FIRST TERM
Dr1 = DT * vkk;

% SECOND TERM
if norm_w_avg < 0.1
    %small omega solution
    Dr2 =  quaternions_to_rot_matrix(qk_1k)'...
        * [ 1/2*DT^2*eye(3)  +  NCP_w_avg  * 1/3*DT^3  +  NCP_w_avg * NCP_w_avg  * 1/8*DT^4 ]...
        * dot_v_avg;
else
    %normal solution
    Dr2 =  quaternions_to_rot_matrix(qk_1k)'...
        * [ 1/2*DT^2*eye(3)  -  NCP_w_avg * 1/norm_w_avg^3*(norm_w_avg*DT*cos(norm_w_avg*DT) - sin(norm_w_avg*DT)) ...
        + NCP_w_avg * NCP_w_avg  * 1/2*1/norm_w_avg^4*( (norm_w_avg*DT)^2 - 2*cos(norm_w_avg*DT) - 2*norm_w_avg*DT*sin(norm_w_avg*DT) + 2 ) ]...
        * dot_v_avg;
end

% THIRD TERM
Dr3 = +1/2*DT^2 * g;

% TOTAL POSITION
rk_1k = rkk + Dr1 + Dr2 + Dr3;

% /// CHANGE ends here

% position state propagation
xek_1k(14:16,1) = rk_1k;

%
%-----------------------------------------------------------------------
%///////////////////////////////////////////////////////////////////////
%-----------------------------------------------------------------------
%
%    C O V A R I A N C E    E S T I M A T E    P R O P A G A T I O N
%
%-----------------------------------------------------------------------
%///////////////////////////////////////////////////////////////////////
%-----------------------------------------------------------------------
%
% Useful Quantities

% /// CHANGE begins here (Nikolas, July 9, 2005) (moved this portion further
% up in the code, to use these quantities also during vel/pos integration)
% % norm of avg. rotational velocity
% norm_w_avg = sqrt(w_avg' * w_avg);
% 
% % negative cross product matrix derived using GYRO measurements
% [NCP_w_avg]=neg_cross_product(w_avg);

% /// CHANGE ends here

% // CHANGE starts here // (Nikolas, July 9 2005) (calculated further up!)
% % Previous UNBIASED acceleration estimate - tk
% dot_v_kk = ACCEL_previous - accel_bias_kk;
% 
% % Current UNBIASED acceleration estimate - tk+1
% dot_v_k_1k = ACCEL_current - accel_bias_k_1k;
% 
% % Average UNBIASED acceleration estimate (w/ ACCELERATION)
% dot_v_avg = 0.5 * (dot_v_kk + dot_v_k_1k);
% // CHANGE ends here //

% negative cross product matrix derived using ACCEL measurements
[NCP_dot_v_avg]=neg_cross_product(dot_v_avg);

% ROTATIONAL MATRIX FROM LOCAL TO GLOBAL - 
% The quaternion estimated by the filter is the one from global to
% local coordinates
% Thus the C_q matrix describes transformations from global to
% local coordinates
% The TRANSPOSE of this matrix: LOCAL -> GLOBAL
[C_q]=quaternions_to_rot_matrix(qkk);

%-----------------------------------------------------------------------
% STATE (SYSTEM) MATRIX CALCULATIONS  --- Fk --- ( 15 x 15 )
%-----------------------------------------------------------------------


% Computation of ZERO SUB-MATRIX Elements
F13 = zeros(3,3);
F14 = zeros(3,3);
F15 = zeros(3,3);

F21 = zeros(3,3);
F23 = zeros(3,3);
F24 = zeros(3,3);
F25 = zeros(3,3);

F35 = zeros(3,3);

F41 = zeros(3,3);
F42 = zeros(3,3);
F43 = zeros(3,3);
F45 = zeros(3,3);

% Computation of IDENTITY SUB-MATRIX Elements
F22 = eye(3);
F33 = eye(3);
F44 = eye(3);
F55 = eye(3);

% Computation of NON-TRIVIAL TERMS (time dependent - GYRO, ACCEL)
% F11 ----------------------------------------------------------
F11 = expm(NCP_w_avg*DT);

% F12 ----------------------------------------------------------
if norm_w_avg < 0.03
   F12 = - [eye(3) * DT + NCP_w_avg * (DT^2)/2 + NCP_w_avg * NCP_w_avg * (DT^3)/6] ;
else
   F12 =  -[ eye(3) * DT + NCP_w_avg * ( 1 - cos(norm_w_avg*DT) )/(norm_w_avg^2) + NCP_w_avg * NCP_w_avg * ( norm_w_avg*DT - sin(norm_w_avg * DT) )/(norm_w_avg^3) ] ;
end

% F31 ----------------------------------------------------------
F31 = - C_q' * NCP_dot_v_avg * F12;

% F32 ----------------------------------------------------------
if norm_w_avg < 0.5
  F32 = - C_q' * NCP_dot_v_avg * [eye(3) * (DT^2)/2 + NCP_w_avg * ...
		    (DT^3)/6 + NCP_w_avg * NCP_w_avg * (DT^4)/24];
else
  F32 = - C_q' * NCP_dot_v_avg * [eye(3) * (DT^2)/2 + NCP_w_avg * ( ...
      norm_w_avg*DT - sin(norm_w_avg*DT) )/(norm_w_avg^3) + NCP_w_avg ...
		    * NCP_w_avg * ( cos(norm_w_avg*DT) - 1 + ...
				    (norm_w_avg^2 * DT^2)/2 )/(norm_w_avg^4)];
end

% F34 ----------------------------------------------------------
F34 = - C_q' * DT;

% F51 ----------------------------------------------------------
F51 = - F32;

% F52 ----------------------------------------------------------
if norm_w_avg < 0.5
  F52 = - C_q' * NCP_dot_v_avg * [eye(3) * (DT^3)/6 + NCP_w_avg * ...
		    (DT^4)/24 + NCP_w_avg * NCP_w_avg * (DT^5)/120];
else
  F52 = - C_q' * NCP_dot_v_avg * [eye(3) * (DT^3)/6 + NCP_w_avg * ( ...
      cos(norm_w_avg*DT) - 1 + (norm_w_avg^2 * DT^2)/2 )/(norm_w_avg^4) + NCP_w_avg * NCP_w_avg * ( sin(norm_w_avg*DT) - norm_w_avg*DT + (norm_w_avg^3 * DT^3)/6 )/(norm_w_avg^5)];
end

% F53 ----------------------------------------------------------
F53 = eye(3) * DT;

% F54 ----------------------------------------------------------
F54 = - 0.5 * C_q' * (DT^2);

% putting all together
Fk = [F11 F12 F13 F14 F15;
      F21 F22 F23 F24 F25;
      F31 F32 F33 F34 F35;
      F41 F42 F43 F44 F45;
      F51 F52 F53 F54 F55];
      
%-----------------------------------------------------------------------
% STATE NOISE COVARIANCE MATRIX CALCULATIONS   --- Qdk --- ( 15 x 15 )
%-----------------------------------------------------------------------

% Computation of ZERO SUB-MATRIX Elements
Qd41=zeros(3,3);
Qd42=zeros(3,3);

Qd14=zeros(3,3);
Qd24=zeros(3,3);

% Computation of NON-TRIVIAL TERMS (time dependent - GYRO, ACCEL)
% QD11 -----------------------------------------------------------------
threshold = 1;
if norm_w_avg<threshold
  Qd11 = sigma_r^2 * DT * eye(3) + sigma_w^2 * [ eye(3) * (DT^3)/3 ...
		    + NCP_w_avg * NCP_w_avg * 2* (DT^5)/120];
else
  Qd11 = sigma_r^2 * DT * eye(3) + sigma_w^2 * [ eye(3) * (DT^3)/3 + NCP_w_avg * NCP_w_avg * ( (norm_w_avg^3 * DT^3)/3 +2* sin(norm_w_avg*DT) -2*norm_w_avg*DT)/(norm_w_avg^5)];
end

% QD12 -----------------------------------------------------------------
if norm_w_avg<threshold
  Qd12 = - sigma_w^2 * [ eye(3) * (DT^2)/2 + NCP_w_avg * (DT^3)/6 + NCP_w_avg * NCP_w_avg * (DT^4)/24 ];
else
  Qd12 = - sigma_w^2 * [ eye(3) * (DT^2)/2 + NCP_w_avg * ( ...
      norm_w_avg*DT-sin(norm_w_avg*DT) )/(norm_w_avg^3) + NCP_w_avg ...
		    * NCP_w_avg * ( cos(norm_w_avg*DT) -1 + (norm_w_avg^2 ...
						  * DT^2)/2 )/(norm_w_avg^4) ];
end

% QD22 -----------------------------------------------------------------
Qd22 = sigma_w^2 * DT * eye(3);

% QD13 -----------------------------------------------------------------
if norm_w_avg<threshold
  Qd13 = - ( sigma_r^2 * [ eye(3)*(DT^2)/2 + NCP_w_avg * (DT^3)/6 + ...
		   NCP_w_avg * NCP_w_avg * (DT^4)/24 ] + sigma_w^2 ...
	     * [ eye(3)*(DT^4)/8 + NCP_w_avg * 2*(DT^5)/120 + NCP_w_avg ...
	       * NCP_w_avg * 5*(DT^6)/720] ) * NCP_dot_v_avg * C_q;
else
   Qd13 = - ( sigma_r^2 * ...
      [ eye(3)*(DT^2)/2 + NCP_w_avg * ( norm_w_avg*DT- sin(norm_w_avg*DT) )/(norm_w_avg^3) + ...
         NCP_w_avg * NCP_w_avg * ( (norm_w_avg^2 * DT^2)/2 + cos(norm_w_avg*DT) -1 )/(norm_w_avg^4) ] + ... 
      sigma_w^2 * ...
         [ eye(3)*(DT^4)/8 + NCP_w_avg * ...
            ( 2*norm_w_avg*DT - 3*sin(norm_w_avg*DT) + norm_w_avg*DT*cos(norm_w_avg*DT) )/(norm_w_avg^5) + ...
         NCP_w_avg * NCP_w_avg * ...
      ( (norm_w_avg^4 * DT^4)/8 - (norm_w_avg^2 * DT^2)/2 + cos(norm_w_avg*DT) + norm_w_avg*DT*sin(norm_w_avg*DT) -1 )/(norm_w_avg^6) ] )...
         * NCP_dot_v_avg * C_q;
end


% QD23 -----------------------------------------------------------------
if norm_w_avg<threshold
  Qd23 = sigma_w^2 * [ eye(3)*(DT^3)/6 - NCP_w_avg *(DT^4)/24 + NCP_w_avg * NCP_w_avg * (DT^5)/120] * NCP_dot_v_avg * C_q;
else
% // CHANGE only this line:   Qd23 = sigma_w^2 * [ eye(3)*(DT^3)/6 - NCP_w_avg * ( (norm_w_avg^2 * DT^2) + cos(norm_w_avg*DT) -1 )/(norm_w_avg^4) ...
% new subterm: 0.5*norm_w_avg^2 * DT^2 //
   Qd23 = sigma_w^2 * [ eye(3)*(DT^3)/6 - NCP_w_avg * ( (0.5*norm_w_avg^2 * DT^2) + cos(norm_w_avg*DT) -1 )/(norm_w_avg^4) ...
         + NCP_w_avg * NCP_w_avg * ( sin(norm_w_avg*DT) + (norm_w_avg^3 * DT^3)/6 - norm_w_avg*DT )/(norm_w_avg^5) ] ...
      * NCP_dot_v_avg * C_q;
end


% QD33 -----------------------------------------------------------------
if norm_w_avg<threshold
   Qd33 = - sigma_r^2 * C_q' * NCP_dot_v_avg * [  eye(3) * (DT^3)/3 + NCP_w_avg * NCP_w_avg * 2*(DT^5)/120 ] * NCP_dot_v_avg * C_q ...
      - sigma_w^2 * C_q' * NCP_dot_v_avg * [ eye(3)*(DT^5)/20 + NCP_w_avg * NCP_w_avg * 10*(DT^7)/5040 ] * NCP_dot_v_avg * C_q ...
      + sigma_g^2 * eye(3) * DT + sigma_v^2 * eye(3) * (DT^3)/3;
else
  Qd33 = - sigma_r^2 * C_q' * NCP_dot_v_avg * [  eye(3) * (DT^3)/3 + NCP_w_avg * NCP_w_avg *( 2*sin(norm_w_avg*DT) -2*norm_w_avg*DT + (norm_w_avg^3 * DT^3)/3 )/(norm_w_avg^5)  ] * NCP_dot_v_avg * C_q ...
      - sigma_w^2 * C_q' * NCP_dot_v_avg * [ eye(3)*(DT^5)/20 + NCP_w_avg * NCP_w_avg * ( (norm_w_avg^5 * DT^5)/20 -2*norm_w_avg*DT - (norm_w_avg^3 * DT^3)/3 + 4*sin(norm_w_avg*DT) - 2*norm_w_avg*DT*cos(norm_w_avg*DT) )/(norm_w_avg^7) ] * NCP_dot_v_avg * C_q ...
      + sigma_g^2 * eye(3) * DT + sigma_v^2 * eye(3) * (DT^3)/3;
end


% QD34 -----------------------------------------------------------------
Qd34 = - sigma_v^2 * C_q' * (DT^2)/2;

% QD44 -----------------------------------------------------------------
Qd44 = DT * sigma_v^2 * eye(3);

% QD15 -----------------------------------------------------------------
if norm_w_avg<threshold
   Qd15 = - sigma_r^2 * [ eye(3) * (DT^3)/6 + NCP_w_avg * 2*(DT^4)/24 + NCP_w_avg * NCP_w_avg * 3*(DT^5)/120] * NCP_dot_v_avg * C_q ...
      - sigma_w^2 * [ eye(3) * (DT^5)/30 + NCP_w_avg * 5*(DT^6)/720 + NCP_w_avg * NCP_w_avg * 11 *(DT^7)/5040 ] * NCP_dot_v_avg * C_q;
else
   Qd15 = - sigma_r^2 * ...
      [ eye(3) * (DT^3)/6 + NCP_w_avg * ( -2*cos(norm_w_avg*DT) - norm_w_avg*DT*sin(norm_w_avg*DT) + 2 )/(norm_w_avg^4) + ...
         NCP_w_avg * NCP_w_avg * ( norm_w_avg*DT + (norm_w_avg^3 * DT^3)/6 -2*sin(norm_w_avg*DT) + norm_w_avg*DT*cos(norm_w_avg*DT) )/(norm_w_avg^5) ] * NCP_dot_v_avg * C_q ...
      - sigma_w^2 * ...
         [ eye(3) * (DT^5)/30 + NCP_w_avg * ( (norm_w_avg^2 * DT^2)/2 -2*cos(norm_w_avg*DT) -2*norm_w_avg*DT*sin(norm_w_avg*DT) + ((norm_w_avg^2 * DT^2)/2)*cos(norm_w_avg*DT) + 2 )/(norm_w_avg^6) + ...
            NCP_w_avg * NCP_w_avg * ( 2*norm_w_avg*DT + (norm_w_avg^5 * DT^5)/30 - (norm_w_avg^3 * DT^3)/6 - 4*sin(norm_w_avg*DT) + 2*norm_w_avg*DT*cos(norm_w_avg*DT) + ((norm_w_avg^2 * DT^2)/2)*sin(norm_w_avg*DT) )/(norm_w_avg^7) ] * NCP_dot_v_avg * C_q;
end


% QD25 -----------------------------------------------------------------
if norm_w_avg<threshold
  Qd25 = sigma_w^2 * [ eye(3) * (DT^4)/24 - NCP_w_avg * (DT^5)/120 + NCP_w_avg * NCP_w_avg * (DT^6)/720 ] * NCP_dot_v_avg * C_q ;
else
   Qd25 = sigma_w^2 * [ eye(3) * (DT^4)/24 ...
         - NCP_w_avg * ( sin(norm_w_avg*DT) - norm_w_avg*DT + (norm_w_avg^3 * DT^3)/6 )/(norm_w_avg^5) ...
         + NCP_w_avg * NCP_w_avg * ( -cos(norm_w_avg*DT) - (norm_w_avg^2 * DT^2)/2 + (norm_w_avg^4 * DT^4)/24 + 1 )/(norm_w_avg^6) ] ...
      * NCP_dot_v_avg * C_q ;
end


% QD35 -----------------------------------------------------------------
if norm_w_avg<threshold
   Qd35 = - sigma_r^2 * C_q' * NCP_dot_v_avg * ...
      [ eye(3) * (DT^4)/8 + NCP_w_avg * 2*(DT^5)/120 + NCP_w_avg * NCP_w_avg * 5*(DT^6)/720 ] * ...
      NCP_dot_v_avg * C_q ...
   - sigma_w^2 * C_q' * NCP_dot_v_avg * ...
      [ eye(3) * (DT^6)/72 + NCP_w_avg * 5*(DT^7)/5040 + NCP_w_avg * NCP_w_avg * 21*(DT^8)/40320 ] * ...
      NCP_dot_v_avg * C_q ...
   + ((DT^2)/2) * sigma_g^2 * eye(3) + ((DT^4)/8) * sigma_v^2 * eye(3) ;
% // CHANGE only the line above   + ((DT^2)/2) * sigma_g^2 * eye(3) - ((DT^4)/8) * sigma_v^2 * eye(3) ;
% new subterm: + ((DT^4)/8) * sigma_v^2 * eye(3) //
else
   Qd35 = - sigma_r^2 * C_q' * NCP_dot_v_avg * ...
      [ eye(3) * (DT^4)/8 + ...
         NCP_w_avg * ( 2*norm_w_avg*DT - 3*sin(norm_w_avg*DT) + norm_w_avg*DT*cos(norm_w_avg*DT) )/(norm_w_avg^5) + ...
         NCP_w_avg * NCP_w_avg * ( (norm_w_avg^4 * DT^4)/8 -(norm_w_avg^2 * DT^2)/2 + cos(norm_w_avg*DT) + norm_w_avg*DT*sin(norm_w_avg*DT) - 1 )/(norm_w_avg^6) ] * ...
      NCP_dot_v_avg * C_q ...
      - sigma_w^2 * C_q' * NCP_dot_v_avg * ...
         [ eye(3) * (DT^6)/72 + ...
            NCP_w_avg * ( 2*norm_w_avg*DT + (norm_w_avg^3 * DT^3)/6 -5*sin(norm_w_avg*DT) + 3*norm_w_avg*DT*cos(norm_w_avg*DT) + 0.5*(norm_w_avg^2 * DT^2)*sin(norm_w_avg*DT) )/(norm_w_avg^7) + ...
            NCP_w_avg * NCP_w_avg * ( (norm_w_avg^6 * DT^6)/72 - (norm_w_avg^4 * DT^4)/8 + cos(norm_w_avg*DT) + norm_w_avg*DT*sin(norm_w_avg*DT) - 0.5*(norm_w_avg^2 * DT^2)*cos(norm_w_avg*DT) -1 )/(norm_w_avg^8) ] * ...
            NCP_dot_v_avg * C_q ...
           + sigma_g^2 * ((DT^2)/2) * eye(3) + sigma_v^2 * ((DT^4)/8) * eye(3) ;
% CHANGE only the line above           + sigma_g^2 * ((DT^2)/2) * eye(3) - sigma_v^2 * ((DT^4)/8) * eye(3) ;
% new subterm: + sigma_v^2 * ((DT^4)/8) * eye(3) //
end


% QD45 -----------------------------------------------------------------
Qd45 = - 0.5 * sigma_v^2 * C_q * (DT^3)/3;

% QD55 -----------------------------------------------------------------
if norm_w_avg<threshold
   Qd55 = - sigma_r^2 * C_q' * NCP_dot_v_avg * [ eye(3) * (DT^5)/20 + NCP_w_avg * NCP_w_avg * 10*(DT^7)/5040 ] * NCP_dot_v_avg * C_q ...
      - sigma_w^2 * C_q' * NCP_dot_v_avg * [ eye(3) * (DT^7)/252 + NCP_w_avg * NCP_w_avg * 42*(DT^9)/362880 ] * NCP_dot_v_avg * C_q ...
      + sigma_g^2 * (DT^3)/3 * eye(3) + sigma_v^2 * (DT^5)/20 * eye(3) ;
else
   Qd55 = - sigma_r^2 * C_q' * NCP_dot_v_avg * ...
      [ eye(3) * (DT^5)/20 + NCP_w_avg * NCP_w_avg * ( (norm_w_avg^5 * DT^5)/20 - (norm_w_avg^3 * DT^3)/3 - 2*norm_w_avg*DT + 4*sin(norm_w_avg*DT) - 2*norm_w_avg*DT*cos(norm_w_avg*DT) )/(norm_w_avg^7) ] * ...
      NCP_dot_v_avg * C_q ...
   - sigma_w^2 * C_q' * NCP_dot_v_avg * ...
      [ eye(3) * (DT^7)/252 + NCP_w_avg * NCP_w_avg * ( (norm_w_avg^7 * DT^7)/252 - (norm_w_avg^5 * DT^5)/20 - 2*norm_w_avg*DT + 6*sin(norm_w_avg*DT) - 4*norm_w_avg*DT*cos(norm_w_avg*DT) - (norm_w_avg^2 * DT^2)*sin(norm_w_avg*DT) )/(norm_w_avg^9) ] * ...
      NCP_dot_v_avg * C_q ...
      + sigma_g^2 * (DT^3)/3 * eye(3) + sigma_v^2 * (DT^5)/20 * eye(3) ;
end


% Computation of TRANSPOSED ELEMENTS
Qd21 = Qd12';
Qd31 = Qd13';
Qd32 = Qd23';
Qd43 = Qd34';
Qd51 = Qd15';
Qd52 = Qd25';
Qd53 = Qd35';
Qd54 = Qd45';

% putting all together
Qdk = [Qd11 Qd12 Qd13 Qd14 Qd15;
       Qd21 Qd22 Qd23 Qd24 Qd25;
       Qd31 Qd32 Qd33 Qd34 Qd35;
       Qd41 Qd42 Qd43 Qd44 Qd45;
       Qd51 Qd52 Qd53 Qd54 Qd55];

%-----------------------------------------------------------------------
% COVARIANCE PROPAGATION EQUATION OF THE KALMAN FILTER
%-----------------------------------------------------------------------

%// CHANGE Tassos 03/14/07
% initialize the matrix
Pk_1k = Pkk;
% put in the appropriate values
Pk_1k(1:15,1:15) = Fk * Pkk(1:15,1:15) * Fk' + Qdk;
% if a copied state exists, then we also need these....
Pk_1k(16:end,1:15) = Pkk(16:end,1:15) * Fk';
Pk_1k(1:15,16:end) = Fk * Pkk(1:15,16:end);



end 

function [A]=quaternions_to_rot_matrix(q)

% quaternions to rotational matrix

A(1,1) = q(1,1)^2 - q(2,1)^2 - q(3,1)^2 + q(4,1)^2  ;
A(1,2) = 2*( q(1,1)*q(2,1) + q(3,1)*q(4,1) ); 
A(1,3) = 2*( q(1,1)*q(3,1) - q(2,1)*q(4,1) );

A(2,1) = 2*( q(1,1)*q(2,1) - q(3,1)*q(4,1) );
A(2,2) = -q(1,1)^2 + q(2,1)^2 - q(3,1)^2 + q(4,1)^2;
A(2,3) = 2*( q(2,1)*q(3,1) + q(1,1)*q(4,1) );

A(3,1) = 2*( q(1,1)*q(3,1) + q(2,1)*q(4,1) );
A(3,2) = 2*( q(2,1)*q(3,1) - q(1,1)*q(4,1) );
A(3,3) = -q(1,1)^2 - q(2,1)^2 + q(3,1)^2 + q(4,1)^2;

end
