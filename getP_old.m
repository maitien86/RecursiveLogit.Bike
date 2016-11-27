%   Computes the probability for each state action pair
%%
function P = getP_old(expV, M)
    %U = sparse(zeros(size(M)));
    n = size(M,1);
    I = find(M);
    [nbnonzero, c] = size(I);
    ind1 = zeros(nbnonzero,1);
    ind2 = zeros(nbnonzero,1);
    s = zeros(nbnonzero,1);
    for i=1:nbnonzero
        [k a] = ind2sub(size(M), I(i));
        ind1(i) = k;
        ind2(i) = a;
        s(i) =  M(k,a) * expV(a)/expV(k);
    end    
    P = sparse(ind1, ind2, s, n, n);
    
end
