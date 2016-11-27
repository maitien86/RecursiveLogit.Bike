%   Initialize optimization structure
%%
function [] = initialize_optimization_structure()
    global Op;
    global isLinkSizeInclusive;
    global isFixedUturn;
    global nbobs;    
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.ETA1 = 0.05;
    Op.ETA2 = 0.75;
    Op.maxIter = 150;
    Op.k = 0;
    Op.n = 14; %number of parameters except link size
    %makes +1 if link size included:
    if isLinkSizeInclusive == true;
        Op.n = Op.n + 1;
    end
    Op.x = -ones(Op.n,1) * 1.5;
    Op.tol = 1e-6;
    Op.radius = 1.0;
    Op.Ak = eye(Op.n);
    Op.H = eye(Op.n);
    Op.Hi = objArray(nbobs);
    for i = 1:nbobs
        Op.Hi(i).value = eye(Op.n);
    end
end
