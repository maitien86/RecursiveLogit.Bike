%Find observations that do forbidden state transitions (incidenceFull=0)
%mostly this occurs when they use a link that had a direction restriction
%FalseObs(i,j)=1 means that on segment j of trip i, there is a forbidden
%state transition from link j to link j+1.

function FalseObs = FindFalseObs()
    file = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/Incidence.txt';
    file_obs = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/TripsObservations.txt';
    incidence = spconvert(load(file));
    Obs = spconvert(load(file_obs));
    FalseObs=zeros(length(Obs(:,1)));
    for i=1:length(Obs(:,1))
        numlinks=find(Obs(i,:));
        for j = 2:(length(numlinks)-1)
            if incidence(Obs(i,j),Obs(i,j+1))==0
                FalseObs(i,j)=1;
            end
        end
    end
end
