% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL_test()

   global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global SampleObs;
    % ----------------------------------------------------
    % If Link size is included
    mu = 1; 
    % MU IS NORMALIZED TO ONE
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    
    if isempty(SampleObs)
        SampleObs = 1:nbobs;
    end
    sample = SampleObs;
    
    Mfull = getM(Op.x);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x);
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
    AttLc = objArray(Op.n);
    for i = 1:Op.n
        AttLc(i).value =  (Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
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
    if minele < -1e-5 %|| minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
    end 
    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE? 
    D = (A * Z - B);
    resNorm = norm(D(:));
    if resNorm > OptimizeConstant.RESIDUAL
       expVokBool = 0;
    end
    if (expVokBool == 0)
        %resNorm
        LL = OptimizeConstant.LL_ERROR_VALUE;
        grad = ones(Op.n,1);
        disp('The parameters not fesible')
        return; 
    end
    % Get gradient
    gradExpV = objArray(Op.n);
    for i = 1:Op.n
        u = M .* (AttLc(i).value); 
        v = sparse(u * Z); 
        p = Atts(i).value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
        p(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
        p = sparse(p);
        gradExpV(i).value =  sparse(A\(v + p)); 
    end
    
    % Compute the LL and gradient.
    
    for n = 1:nbobs    
%        n
        dest = Obs(sample(n), 1);
        orig = Obs(sample(n), 2);
        expV = Z(:,dest - lastIndexNetworkState);
        expV = full(abs(expV));         
        lnPn = - 1 * (1/mu) * log(expV(orig));
        for i = 1: Op.n
            Gradient(n,i) = - gradExpV(i).value(orig,dest - lastIndexNetworkState)/expV(orig);
        end
        
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        
        path = Obs(sample(n),:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1));
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + Atts(j).value(path(i),path(i+1));
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

%%
