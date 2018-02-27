function B = AccumulateCompactToFull(A, Ac, idx_row, idx_col)

B = A;
block_row_size = size(Ac,1) / length(idx_row);
block_col_size = size(Ac,2) / length(idx_col);

if size(idx_row,2) == 1
    idx_row = idx_row';
end

for k = 1:length(idx_row)
    for l = 1:length(idx_col)
        B(block_row_size*(idx_row(k)-1)+1:block_row_size*idx_row(k), ...
          block_col_size*(idx_col(l)-1)+1:block_col_size*idx_col(l)) = ...
        A(block_row_size*(idx_row(k)-1)+1:block_row_size*idx_row(k), ...
          block_col_size*(idx_col(l)-1)+1:block_col_size*idx_col(l)) + ...
        Ac(block_row_size*(k-1)+1:block_row_size*k, ...
           block_col_size*(l-1)+1:block_col_size*l);
    end
end