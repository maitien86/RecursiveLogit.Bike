%   Generate obs
%%
global Op;
global Obs;

disp('Observation generating ....')
Op.n = 5;
x0 = [-2.0 -1.0 -1.0 -20.0 -5 ]';
ODpairs = Obs(:,1:2);
nbobsOD = 1;
filename = './simulatedData/ObsGLS.txt';
generateObs(filename, x0, ODpairs, nbobsOD);
