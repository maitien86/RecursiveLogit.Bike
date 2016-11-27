function [] = Estimation_IMtest(nVariable, fileObs) 

    globalVar;
    global isLinkSizeInclusive;
    global file_observations;
    global isFixedUturn;
    tic;
    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    
    if str2num(nVariable) == 4
        isLinkSizeInclusive = false;
    else
        isLinkSizeInclusive = true;
    end
    isFixedUturn = true;
    file_observations = fileObs;
    
    loadData;
    
    Op = Op_structure;
    initialize_optimization_structure();
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BHHH;
    Gradient = zeros(nbobs,Op.n);
    
    if isLinkSizeInclusive == false
        Op.x = [-1.8, -1.0, -1.0]';
    else
        Op.x = [-1.8, -1.0, -1.0, -1.20]';
    end

    % Optimize with full choice set
    EstimationScript;
    
    TXT = ['RL|',fileObs,'|', num2str(Op.n),'|'];
    Op.x
    %TXT = [TXT,IMtest(Op.x)];
    TXT = [TXT,IM_fulltest(Op.x),'\n'];
    fileID = fopen('IMtestResults.txt','at+');
    fprintf(fileID,TXT);
    fclose(fileID);
    
    fileID = fopen('ComputationalTime.txt','at+');
    fprintf(fileID,['IM|',TXT,'Elapsed time:', num2str(toc),'\n']);
    fclose(fileID);
end
