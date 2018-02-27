function [idx] = SUTrilInd(N,i_,j_)
    i = min(i_, j_);
    j = max(i_, j_);
    idx = 1/2*(i-1)*(2*N-i)+ j-i;
end