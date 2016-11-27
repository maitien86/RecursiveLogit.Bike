%% Data processing
%% Reading the data
dataLinkIDs = xlsread('./Input/LinkAttributesPortland.xlsx');

%% Link incidence
[nLinks,nbAtts] = size(dataLinkIDs);
nbiLinks = nLinks;

%Zone to Node
filename='./Input/ZoneToNode.xlsx';
Table=readtable(filename);
CentroidToNode=table2array(Table(:,2:3));
NbOD=length(CentroidToNode(:,1));

% Match each OD node to links
DestLinks=cell(1,NbOD);
OrigLinks=cell(1,NbOD);
for k=1:NbOD
    DestLinks{k}=find(dataLinkIDs(:,3)==CentroidToNode(k,2));
    OrigLinks{k}=find(dataLinkIDs(:,2)==CentroidToNode(k,2));
end

%% Creating incidence matrix ###
nStartLinks = double(NbOD);
nAbsorbingLinks = double(NbOD);
A = dataLinkIDs(:,3); %To Nodes (end of link)
B = dataLinkIDs(:,2); %From Nodes (start of link)
Incidence = sparse(nbiLinks+nStartLinks,nbiLinks+nStartLinks+nAbsorbingLinks)
for i = 1:nbiLinks
    U = find(B == A(i));%for each link i, look at all the links U that connect to the End Node.
    Incidence(i,U) = 1;% this means we can go from link i to link U
end
for i = 1: nAbsorbingLinks
    Incidence(DestLinks{i},nbiLinks + nStartLinks + i) = 1;
end
for i = 1: nStartLinks
    Incidence(nbiLinks + i,OrigLinks{i}) = 1;
end

[i,j,val] = find(Incidence);
data_dump = [i,j,val];
save('IncidencePortlandNetwork.txt','data_dump','-ascii');

%dlmwrite('fullLinksAttsPortlandNetwork.txt',dataFullLinksIDs);