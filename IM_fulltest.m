%   Information matrix test (White - 1982)
%%
function testResults = IM_fulltest(beta)
    globalVar;    
    tic;
    step = 1e-9;
    Op.x = beta;
    vsize = Op.n * (Op.n + 1) / 2;
    [~,~,Hes,Hs] = getAnalyticalHessian();
    %LL(beta);    
    h = eye(Op.n) * step;
    d = zeros(nbobs , vsize);
    for n = 1: nbobs
        h1 = Gradient(n,:)' * Gradient(n,:);
        h2 = Hs(n).value;        
        d(n,:) = h1(triu(true(size(h1))))' - h2(triu(true(size(h2))))';        
    end
    H_st = BHHH();
    deltaD = zeros(vsize, Op.n);
    grad = Gradient;
    invA = inv(Hes);
    Hdiff = H_st - Hes;
    D = Hdiff(triu(true(size(Hdiff))));
    Op.x = beta;
    for j= 1 : Op.n
        Op.x = (beta + h(:,j));
        [~,~,A1,~] = getAnalyticalHessian();
        Hdiff = BHHH() - A1; 
        D1 = Hdiff(triu(true(size(Hdiff))));
        deltaD(:,j) = (D1-D)/step;
    end
    V = zeros(vsize);
    U = objArray(Op.n);
    for i = 1: Op.n
        U(i).value = 0;
    end
    for n = 1: nbobs
        v = d(n,:)' + deltaD * invA * (grad(n,:))';
        V =  V + (v * v' - V)/n;
        for i = 1:Op.n
            u = d(n,i)' + deltaD(i,:) * invA * (grad(n,:))';
            U(i).value =  U(i).value + (u * u' - U(i).value)/n;
        end
    end
    testValue = zeros(Op.n + 1,1);
    V
    testValue(1) = nbobs * D' * inv(V) * D;
    for i = 1: Op.n
        testValue(i+1) =  nbobs * D(i,:)' * inv(U(i).value) * D(i,:);
    end
    testValue(1);
    testResults = sprintf('RL|%10.5e|%7d|', testValue(1), vsize); 
    
    %fprintf(' Chi test value : %e \n', testValue(1));
    %fprintf(' df             : %d \n', vsize);
    %fprintf(' p - value      : %e \n', 1 - chi2cdf(testValue(1),vsize));
    %fprintf(' Elapsed time   : %e (seconds) \n', toc);
end