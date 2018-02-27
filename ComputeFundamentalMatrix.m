function F = ComputeFundamentalMatrix(uvmatrix)

    U = uvmatrix(:,1:2);
    V = uvmatrix(:,3:4);
    
    A = [];
    for i = 1:size(uvmatrix,1)
        A = [A; reshape([U(i,:) 1]'*[V(i,:) 1], [1,9])];
    end
        
    if(rank(A) < 8)
        F = [];
        return;
    end
    [~, ~, V] = svd(A);
    F = reshape(V(:,end), [3,3])';
    
    [Uf,Sf,Vf] = svd(F);
    if(Sf(2,2)/Sf(1,1) < 5e-7)
        F = [];
        return
    end
    Sf(3,3) = 0;
    
    F = Uf*Sf*Vf';
end