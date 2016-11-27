function [] = Hausman_IIA_test(nVariable, fileObs) 

    globalVar;
    global isLinkSizeInclusive;
    global file_observations;
    global isFixedUturn;
    global incidenceFull;
    global Obs;
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
%     % Optimize with full choice set
    EstimationScript;
    % Get Cov matrix
    x_full = Op.x;
    fprintf('\n Getting analytical Hessian ....\n');
    Hessian = getHessian();
    InvH = inv(Hessian);
    CoV_full = (InvH * BHHH() * InvH)/nbobs;
    
    %% Optimizing with sub choice set
    %nbobs = 1000;
    %[incidenceFull, Obs] = CreatSubNetWork(500, incidenceFull, Obs);
    incidenceFull = spconvert(load('./Input/linkIncidence_sub1000.txt'));
    % Resample Obs
    i = 1;
    while(i <= size(Obs,1))
        path = Obs(i,:);
        if isPathExiste(path) == false
            Obs(i,:) = [];
            i = i-1;
        end
        i = i+1;
    end    
    nbobs = size(Obs,1);
    Atts = getAtt();
    Op.radius = 1.0;
    Gradient = zeros(nbobs,Op.n);
    EstimationScript;
    x_sub = Op.x;
    fprintf('\n Getting analytical Hessian ....\n');
    Hessian = getHessian();
    InvH = inv(Hessian);
    CoV_sub = (InvH * BHHH() * InvH)/nbobs;
    
    diffCov = -(CoV_full - CoV_sub)
    testVl = (x_full - x_sub)' * inv(diffCov) * (x_full - x_sub)
    TXT = sprintf('%10.5e|%7d|\n', testVl, Op.n);    
    TXT = ['RL|',fileObs,'|', num2str(Op.n),'|', TXT];
    fileID = fopen('IIAtestResults.txt','at+');
    fprintf(fileID,TXT);
    fclose(fileID);

    fileID = fopen('ComputationalTime.txt','at+');
    fprintf(fileID,['IIA|',TXT,'Elapsed time:', num2str(toc),'\n']);
    fclose(fileID);
end

function [ok] = isPathExiste(path)
    global incidenceFull;
    ok = true;
    for j=2: size(find(path),2)-1
        if incidenceFull(path(j), path(j+1)) == 0
            ok = false;
            return;
        end
    end

end
