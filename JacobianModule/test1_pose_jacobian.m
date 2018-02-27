% Data generation
for trials = 1:10000
close all;
clear all;

% orientation:
Ci_q_Cr = randn(4,1); Ci_q_Cr = Ci_q_Cr/norm(Ci_q_Cr);
if Ci_q_Cr(4) < 0
    Ci_q_Cr(4) = -Ci_q_Cr(4);
end
Ci_R_Cr = quat2rot(Ci_q_Cr);
% position:
Cr_p_Ci = 10*randn(3,1);

% Model:
% z > 0!!
featureNums = 20;
Ci_p_f = 5*randn(3, featureNums);
Ci_p_f(:, Ci_p_f(3,:)<0) = -Ci_p_f(:, Ci_p_f(3,:)<0);
Ci_p_f(3,:) = Ci_p_f(3,:) + 5; % good condition

% Geometric model
Cr_p_f = Ci_R_Cr'*Ci_p_f + repmat(Cr_p_Ci, [1 featureNums]);

% Projection onto image
Ci_z_f = Ci_p_f(1:2,:) ./ repmat(Ci_p_f(3,:), [2,1]);

% Pertubation the state vector Cr_p_Ci and Ci_R_Cr
std_p = 1.5;
std_tt = 0.1 ;
std_f = 0.3;

tt_pert = std_tt*randn(3,1);
dq = [1/2*tt_pert; sqrt(1 - 1/4*norm(tt_pert)^2)];
Ci_q_Cr_hat = quat_mul(dq, Ci_q_Cr);
Ci_R_Cr_hat = quat2rot(Ci_q_Cr_hat);
Cr_p_Ci_hat = Cr_p_Ci + std_p*randn(3,1);

% Testing Jacobian
Ci_p_f_hat = Ci_R_Cr_hat*(Cr_p_f - repmat(Cr_p_Ci_hat, [1 featureNums]));
Ci_z_f_hat = Ci_p_f_hat(1:2,:) ./ repmat(Ci_p_f_hat(3,:), [2,1]);

[~, res] = dprojection_dp_dtt(Cr_p_f, Ci_R_Cr_hat, Cr_p_Ci_hat, Ci_z_f);
[Jpose, ~] = dprojection_dp_dtt(Cr_p_f, Ci_R_Cr, Cr_p_Ci, Ci_z_f);

norm(res)
norm(-res - Jpose*[Cr_p_Ci_hat-Cr_p_Ci;tt_pert])


% Gauss-Newton iteration
Cr_p_Ci_iter = Cr_p_Ci_hat;
Ci_q_Cr_iter = Ci_q_Cr_hat;
iter = 1;
norm_res_log = [norm(res)];
while norm(res) > 1e-5
    [Jpose, res] = dprojection_dp_dtt(Cr_p_f, quat2rot(Ci_q_Cr_iter), Cr_p_Ci_iter, Ci_z_f);
    
    % Solver:
    dx = (Jpose'*Jpose) \ (Jpose'*res);
    dp = dx(1:3); dtt = dx(4:6);
    
    % update:
    Cr_p_Ci_iter = Cr_p_Ci_iter + dp;
    dq = [1/2*dtt; sqrt(1 - 1/4*norm(dtt)^2)];
    Ci_q_Cr_iter = quat_mul(dq, Ci_q_Cr_iter);
    
    % Logger
    norm_res_log = [norm_res_log norm(res)];
end

[Cr_p_Ci_hat Cr_p_Ci_iter Cr_p_Ci]
[Ci_q_Cr_hat Ci_q_Cr_iter Ci_q_Cr]

figure;
semilogy(norm_res_log);
pause;
end