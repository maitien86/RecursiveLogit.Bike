function [] = prediction(instanceN)
instanceNum=str2double(instanceN);
%prediction predicts LL loss for subsample of the data
%   Detailed explanation goes here

Credits;
% Initialize email notification
globalVar;

notifyMail('set','maelle.zimmermann@gmail.com','sntal2908');
global resultsTXT; 
global AttsNames;
global isLinkSizeInclusive;
global LSatt; %load different LSatt depending on network!
global nbobs; %before this was ommitted from code??
isLinkSizeInclusive = true;

%change for loadData_2WaysNetwork if we want the other input data (bi links
%everywhere in the network)
%loadData_2WaysNetwork;
loadDataNRL;
Op = Op_structure;
initialize_optimization_structure();

Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
Op.Hessian_approx = OptimizeConstant.BFGS;
Gradient = zeros(nbobs,Op.n);

% Prediction parameters
global SampleObs
PredSample = spconvert(load('/home/zimmae/Documents/RL2.6.DeC.GLS/PredSample648Obs.txt'));
SampleObs = PredSample(instanceNum*2-1,:);
nbobs = size(find(SampleObs),2);

%---------------------------
%Starting optimization
tic ;
%progTest
disp('Start Optimizing ....')

% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','GradObj','on');
% [x,fval,exitflag,output,grad] = fminunc(@LL,Op.x,options)

if isLinkSizeInclusive == true
    [Op.value, Op.grad ] = getLL(); %write getLL() for link size and getLL_test() else
else
    [Op.value, Op.grad ] = getLL_test();
end
PrintOut(Op);
% print result to string text
header = [sprintf('%s \n',file_observations) Op.Optim_Method];
header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
resultsTXT = header;
%------------------------------------------------
while (true)    
  Op.k = Op.k + 1;
  if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
    ok = line_search_iterate();
    if ok == true
        PrintOut(Op);
    else
        disp(' Unsuccessful process ...')
        break;
    end
  else
    ok = btr_interate();
    PrintOut(Op);
  end
  [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
  %----------------------------------------
  % Check stopping criteria
  if(isStop == true)
      isSuccess
      fprintf('The algorithm stops, due to %s', Stoppingtype);
      resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
      break;
  end

end

% %Save to file
% formatout1='yymmdd';
% formatout2='hhMM';
% date1=datestr(now,formatout1);
% date2=datestr(now,formatout2);
% FolderString=horzcat('./Results/',date1);
% if ~exist(FolderString, 'dir')
%   % The folder does not exist.
%   % Create that folder.
%   mkdir(FolderString);
% end
% FileString=horzcat('./Results/',date1,'/',date2,'.txt');
% FileID=fopen(FileString,'w');
% fprintf(FileID,'The algorithm stops, due to %s \n', Stoppingtype);
% fprintf(FileID,'The attributes are \n');
% for cellitem = 1:length(AttsNames)
%     fprintf(FileID,'%s \n',AttsNames{cellitem});
% end           
% fprintf(FileID,'[Iteration]: %d\n', Op.k);
% fprintf(FileID,'     LL = %f\n', Op.value);
% fprintf(FileID,'     x = \n');
% fprintf(FileID,'         %i\n', Op.x');
% fprintf(FileID,'     norm of step = %f\n', norm(Op.step));
% fprintf(FileID,'     radius = %f\n', Op.delta);  
% fprintf(FileID,'     Norm of grad = %f\n', norm(Op.grad));
% relatice_grad = relative_gradient(Op.value, Op.x, Op.grad, 1.0);
% fprintf(FileID,'     Norm of relative gradient = %f\n', relatice_grad);
% fprintf(FileID,'     Number of function evaluation = %f\n', Op.nFev);
% 
% %Finishing ...
% ElapsedTime = toc
% resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
% resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTime)];
% 
% %Continue save to file
% fprintf(FileID,'\n Number of function evaluation %d \n', Op.nFev);
% fprintf(FileID,'\n Estimated time %d \n', ElapsedTime);
% fprintf(FileID,'Estimated : %5.8f \n',Op.x);
% fprintf(FileID,'Standard deviation : %5.8f \n',Stdev);
% fprintf(FileID,resultsTXT);
% fclose(FileID);

% Predicts LL loss
TXT = 'Prediction:';
TXT = [TXT sprintf('\n Link size = %d \n', isLinkSizeInclusive)];
SampleObs=PredSample(instanceNum*2,:);
nbobs = size(find(SampleObs),2);
if isLinkSizeInclusive == true
    [LL,grad] = getLL(); %Write getLL_test() without LS!
else 
    [LL,grad] = getLL_test();
end
TXT = [TXT sprintf('%d : %f \n', instanceNum, LL)];

%Save to file
    formatout1='yymmdd';
    formatout2='hhMM';
    date1=datestr(now,formatout1);
    date2=datestr(now,formatout2);
    instance=num2str(instanceNum);
    FolderString=horzcat('/home/zimmae/Documents/RL2.6.DeC.GLS/Preds/',date1);
    if ~exist(FolderString, 'dir')
    % The folder does not exist.
    % Create that folder.
    mkdir(FolderString);
    end
    FileString=horzcat('/home/zimmae/Documents/RL2.6.DeC.GLS/Preds/',date1,'/',date2,'inst',instance,'.txt');
    FileID=fopen(FileString,'w');
    fprintf(FileID,TXT);
    fprintf(FileID,'The estimation results on the training sample were \n');
    fprintf(FileID,'The algorithm stops, due to %s', Stoppingtype);
    fprintf(FileID,'     LL = %f\n', Op.value);
    fprintf(FileID,'     x = \n');
    fprintf(FileID,'         %i\n', Op.x');
    fclose(FileID);

try   
    notifyMail('send', resultsTXT);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end

end

