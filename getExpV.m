% Get e^V(k)
%%
function [expV, boolV0] = getExpV(M)
    boolV0 = 1;
    [n m]= size(M);
    b = zeros(n,1);
    b(n)=1;
    b_sp = sparse(b);    
    I = speye(size(M));  
    A = I - M; 
    Z = A\b_sp;
    % Z = linsolve(A,b_sp)   
%     Z(Z == 0) = realmin;
    
    %test "values smaller than 0"
    minele = min(Z);
    if minele == 0 || minele < OptimizeConstant.NUM_ERROR
        boolV0 = 0;
    end  
    
    %set values smaller than real min to realmin
    Z(Z<realmin) = realmin;
    
    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE?  
    
    %test "still solution of system after making neg elements positive"
    resNorm = norm(A * Zabs - b_sp);
    if resNorm > OptimizeConstant.RESIDUAL
        boolV0 = 0;
    end
    
%     %test "is fixed point"
%     Zprime = M*Zabs + b_sp;
%     Zdiff=Zprime-Zabs;
%     AbsZdiff=abs(Zdiff);
%     Diff=AbsZdiff./Zabs;
%     Diff2=AbsZdiff./Zprime;
%     Zbis=Zprime;
%     Zbis(find(Zprime<realmin))=realmin;

    expV = full(Zabs);

end
