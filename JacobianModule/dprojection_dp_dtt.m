function [Jpose,residual] = dprojection_dp_dtt(Cr_p_f, Ci_R_Cr, Cr_p_Ci, Ci_z_f)

% dCr_p_Ci | dCi_tt_Cr

featureNums = size(Cr_p_f,2);
Ci_p_f = Ci_R_Cr*(Cr_p_f - repmat(Cr_p_Ci, [1 featureNums]));
Jpose = zeros(2*featureNums, 6);

Ci_z_f_hat = Ci_p_f(1:2,:) ./ repmat(Ci_p_f(3,:), [2,1]);
residual = reshape(Ci_z_f, [2*featureNums 1]) - reshape(Ci_z_f_hat, [2*featureNums 1]);

for k = 0:featureNums-1
    Jpose(2*k+1:2*k+2,1:3) = 1/Ci_p_f(3,k+1)*[1 0 -Ci_p_f(1,k+1)/Ci_p_f(3,k+1); ...
                                            0 1 -Ci_p_f(2,k+1)/Ci_p_f(3,k+1)]*(-Ci_R_Cr);
    Jpose(2*k+1:2*k+2,4:6) = 1/Ci_p_f(3,k+1)*[1 0 -Ci_p_f(1,k+1)/Ci_p_f(3,k+1); ...
                                            0 1 -Ci_p_f(2,k+1)/Ci_p_f(3,k+1)]*skewsymm(Ci_R_Cr*(Cr_p_f(:,k+1) - Cr_p_Ci));
end