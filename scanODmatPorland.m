%Read demand files
%Commute=spconvert(load('/home/zimmae/Documents/RL2.6.DeC.GLS/Input/2010_dailyBikeTripsCommute.csv'));
%NonCommute=spconvert(load('/home/zimmae/Documents/RL2.6.DeC.GLS/Input/2010_dailyBikeTripsNonComm.csv'));

Commute=dlmread('/home/zimmae/Documents/RL2.6.DeC.GLS/Input/2010_dailyBikeTripsCommute.csv');
NonCommute=dlmread('/home/zimmae/Documents/RL2.6.DeC.GLS/Input/2010_dailyBikeTripsNonComm.csv');

%load incidence matrix
file_linkIncidence = './Input/IncidencePortlandNetwork.txt';
incidenceFull = spconvert(load(file_linkIncidence));
[nbStatesRealStart,nbStatesRealStartAbsorb] = size(incidenceFull); 
nStartLinks=max(Commute(:,1));
nbiLinks=308520;

%% Write OD demand in terms of link IDs
%ODmatCommute=sparse(size(incidenceFull));
%ODmatNonCommute=sparse(size(incidenceFull));
%[row, col]=find(Commute>0);
%for k=1:length(row)
%    i=row(k);%origin Centroid ID
%    j=col(k);%destination Centroid ID
%    ODmatCommute(nbiLinks+i,nbiLinks+nStartLinks+j)=Commute(i,j);
%end
%[row, col]=find(NonCommute>0);
%for k=1:length(row)
%    i=row(k);%origin Centroid ID
%    j=col(k);%destination Centroid ID
%    ODmatNonCommute(nbiLinks+i,nbiLinks+nStartLinks+j)=NonCommute(i,j);
%end
%ODmat=ODmatCommute+ODmatNonCommute;

nbiLinksVect=nbiLinks*ones(length(Commute(:,1)),1);
nStartLinksVect=nStartLinks*ones(length(Commute(:,1)),1);
%[i,j,val] = find(Incidence);
data_dump = [Commute(:,1)+nbiLinksVect,Commute(:,2)+nbiLinksVect+nStartLinksVect,Commute(:,3)];
save('ODmatCommute.txt','data_dump','-ascii');

%[i,j,val] = find(Incidence);
data_dump = [NonCommute(:,1)+nbiLinksVect,NonCommute(:,2)+nbiLinksVect+nStartLinksVect,NonCommute(:,3)];
save('ODmatNonCommute.txt','data_dump','-ascii');