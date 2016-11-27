%   Information matrix test (White - 1982)
%%
function testResults = IMtest(beta)
    globalVar; 
    tic ;
    step = 1e-9;
    Op.x = beta;
    [f,g,Hes,Hs] = getAnalyticalHessian();
    Hes;
    %LL(beta);    
    h = eye(Op.n) * step;
    d = zeros(nbobs,Op.n);
    for n = 1: nbobs
        d(n,:) = (diag(Gradient(n,:)' * Gradient(n,:),0) - diag(Hs(n).value,0))';        
    end
    H_st = BHHH();
    D1 = zeros(1,Op.n);
    deltaD = zeros(Op.n);
    grad = Gradient;
    invA = inv(Hes);
    D = diag(H_st - Hes,0);
    Op.x = beta;
    for j= 1 : Op.n
        Op.x = (beta + h(:,j));
        [f,g,A1,Hs] = getAnalyticalHessian();
        D1 = diag(BHHH() - A1,0);
        deltaD(:,j) = (D1-D)/step;
    end
    V = zeros(Op.n);
    for n = 1: nbobs
        v = d(n,:)' + deltaD * invA * (grad(n,:))';
        V =  V + (v * v' - V)/n;
    end
    testValue = nbobs * D' * inv(V) * D;
    testResults = sprintf('%10.5e|%7d|%7.3e|', testValue, Op.n, 1 - chi2cdf(testValue,Op.n)); 
    
%    fprintf(' IM diagonal test: \n');
%    fprintf(' Chi test value : %e \n', testValue);
%    fprintf(' p - value      : %e \n', 1 - chi2cdf(testValue,Op.n));
%    fprintf(' Elapsed time   : %e (seconds) \n', toc);
end