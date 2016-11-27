%% Generate training sample and holdout sample based on observations
global nbobs
n  = nbobs;
tn = floor(0.8 * nbobs);
nloops = 20;
PredSample = zeros(2*nloops, tn);
for i = 1 : nloops
    per = randperm(n);
    PredSample(2*i-1,1:tn) = per(1:tn);
    PredSample(2*i,1:n - tn) = per(tn + 1: n);
end
PredSample = sparse(PredSample);
[i,j,val] = find(PredSample);
data_dump = [i,j,val];
save('PredSample648Obs.txt','data_dump','-ascii');

% %% Generate tranning sample and hodout sample based on ODs
% %%Read real Observations file
% file_observations = './Input/observationsForEstimBAI.txt';
% Obs = spconvert(load(file_observations));
% [nbobs, maxstates] = size(Obs);
% A = Obs(:,1:2);
% [C,ia,ic] = unique(A,'rows');
% n = size(ia,1);
% tn = round(n * 0.8);
% PredSample = zeros(40, tn);
% for i = 1 : 20
%     per = randperm(n);
%     train = per(1:tn);
%     test = per(tn + 1: n);
%     trainIdx = (find(ismember(ic,train)));
%     testIndx = (find(ismember(ic,test)));    
%     PredSample(2*i-1,1:size(trainIdx,1)) = trainIdx;
%     PredSample(2*i,1:size(testIndx,1)) = testIndx;
% end
% %PredSample(1,1:2) = 1:2;
% %PredSample(2,1:2) = 3:4;
% PredSample = sparse(PredSample);
% [i,j,val] = find(PredSample);
% data_dump = [i,j,val];
% save('../PredSampleODs.21.40.txt','data_dump','-ascii');