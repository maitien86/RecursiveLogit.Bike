%   Get analytical Hessian matrix
%   Link size is included
%%
function [LL, grad, Hessian, Hs] = getAnalyticalHessian_test()
    
    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global isLinkSizeInclusive;
    global lastIndexNetworkState;
    % ----------------------------------------------------
    % If Link size is included
    mu = 1; 
    % MU IS NORMALIZED TO ONE
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    
    Mfull = getM(Op.x, isLinkSizeInclusive);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x, isLinkSizeInclusive);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    % Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    % Initialize
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end

    % b matrix:
    N = size(M,1);
    b = sparse(zeros(N,1));
    b(N) = 1;
    B = sparse(zeros(N, maxDest - lastIndexNetworkState));
    B(N,:) = ones(1,maxDest - lastIndexNetworkState);
    for i = 1: maxDest - lastIndexNetworkState
        B(1:lastIndexNetworkState,i) = Mfull(:, i+lastIndexNetworkState);
    end
    A = speye(size(M)) - M;
    Z = A\B;
    % Check feasible
    minele = min(Z(:));
    expVokBool = 1;
    if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
    end 
    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE? 
    D = (A * Z - B);
    resNorm = norm(D(:));
    if resNorm > OptimizeConstant.RESIDUAL
       expVokBool = 0;
    end
    if (expVokBool == 0)
            LL = OptimizeConstant.LL_ERROR_VALUE;
            grad = ones(Op.n,1);
            disp('The parameters not fesible')
            return; 
    end
    % Get gradient
    gradExpV = objArray(Op.n);
    gradM = objArray(Op.n);
    for i = 1:Op.n
        u = M .* (AttLc(i).Value); 
        gradM(i).value = u;
        v = sparse(u * Z); 
        p = Atts(i).Value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
        p(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
        p = sparse(p);
        gradExpV(i).value =  sparse(A\(v + p)); 
    end
    % Evaluate hessian
    hessianM = objMatrix(Op.n, Op.n);
    hessianB = objMatrix(Op.n, Op.n);
    hessianZ = objMatrix(Op.n, Op.n);
    Hessian = zeros(Op.n);
    H = zeros(Op.n);
    Hs = objArray(nbobs);
    % Get second order derivative of M
    for i = 1: Op.n
        for j = 1: Op.n
            u = AttLc(i).Value .* AttLc(j).Value;
            v = M .* u;
            hessianM(i,j).value = v;
        end
    end
    % Get second derivative of B   
    for i = 1: Op.n
        for j = 1: Op.n
            p = Atts(i).Value(:,lastIndexNetworkState+1 : maxDest) .*  Atts(j).Value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
            p(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
            hessianB(i,j).value = sparse(p);        
        end
    end
    for i = 1: Op.n
        for j = 1: Op.n
            p1 = hessianM(i,j).value * Z;
            p2 = gradM(i).value * gradExpV(j).value;
            p3 = gradM(j).value * gradExpV(i).value;
            p4 = hessianB(i,j).value;
            hessianZ(i,j).value = sparse(A\(p1 + p2 + p3 + p4));
        end
    end

    % Compute the LL, grad abd hessian
    
    for n = 1:nbobs    
        n
        dest = Obs(n, 1);
        orig = Obs(n, 2);
        expV = Z(:,dest - lastIndexNetworkState);
        expV = full(abs(expV));         
        lnPn = - 1 * (1/mu) * log(expV(orig));
        % gradient
        for i = 1: Op.n
            Gradient(n,i) = - gradExpV(i).value(orig,dest - lastIndexNetworkState)/expV(orig);
        end
        % hessian
        for i = 1: Op.n
            for j = 1: Op.n
                H(i,j) = (hessianZ(i,j).value(orig,dest - lastIndexNetworkState)/expV(orig) - Gradient(n,i) * Gradient(n,j));
            end
        end
        Hs(n).value = H;
        Hessian = Hessian + (H - Hessian)/n;       
        
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        
        path = Obs(n,:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1));
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + Atts(j).Value(path(i),path(i+1));
            end
        end
        Gradient(n,:) = Gradient(n,:) + sumInstX;
        lnPn = lnPn + (1/mu)*sumInstU ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

function [gradExpV] = getGradExpV(M, Att, op, expV)
    I = speye(size(M));  
    A = I - M; 
    for i = 1:op.n
        u = M .* (Att(i).Value); 
        v = u * expV; 
        gradExpV(i) = Matrix2D(A\v); 
    end
end
