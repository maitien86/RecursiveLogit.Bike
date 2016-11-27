% Transfer matrix A to square matrix, add dummies node
function M = transferMarix(A)
   hsize = max(size(A));
   lsize = min(size(A));
   M = sparse(zeros(hsize + 2));
   M(2:size(A,1)+1, 2:size(A,2)+1) = A;
end