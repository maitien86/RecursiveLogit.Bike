global Incidence

%% Get incidence matrix and link attributes
file_linkIncidence = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/Incidence_2WaysNetwork.txt'; 
file_linkAttributes = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/fullLinksAtts_2WaysNetwork.txt';
incidenceFull = spconvert(load(file_linkIncidence));
linkAttributes = load(file_linkAttributes) ; %each link results in two states! (two ways)
[nbRealStates, nbTotalStates] = size(incidenceFull); 

%% Read demand matrix
filename='./Input/bucket_rounded_bike_new.txt';
fid=fopen(filename);
[Nout]=textscan(fid,'%d%d:%d%d:%d%d:%d%d:%d%d:%d','CollectOutput',1,'EmptyValue',0);
fclose(fid);
out=[Nout{:}];
NbOD=out(end,1);
ODdemand=zeros(NbOD);
%Shape into demand matrix indexed with with centroid ID
for i=1:length(out(:,1))
    for j=1:5
        if out(i,2*j) ~= 0
            ODdemand(out(i,1),out(i,2*j))=out(1,2*j+1);
        end
    end
end

%% Translate CentroidID to NodeID
filename='./Input/TAZ_Node_Link.xlsx';
Table=readtable(filename);
CentroidToNode=table2array(Table(:,1:2));
CentroidToNode(end+1,:)=[24,0]; %row 24 missing
CentroidToNode2=sortrows(CentroidToNode);
NbOD=length(CentroidToNode(:,1));

% Match each OD node to links
DestLinks=cell(1,NbOD);
OrigLinks=cell(1,NbOD);
for k=1:NbOD
    DestLinks{k}=find(linkAttributes(:,3)==CentroidToNode2(k,2));
    OrigLinks{k}=find(linkAttributes(:,2)==CentroidToNode2(k,2));
end

%% Data processing: create incidence matrix with appropriate absorbing links
%Reading the data
dataLinkIDs = xlsread('./Input/LinkIDs.xlsx');
[nLinks,nbAtts] = size(dataLinkIDs);
temp = dataLinkIDs(:,4);
temp(find(temp == 2)) = 1;
temp = 2 - temp;

nbiLinks = 2*nLinks;

mappingLinks = zeros(nbiLinks,1); % For mapping
dataFullLinksIDs = zeros(nbiLinks,nbAtts);
dataFullLinksIDs(1:nLinks,:) = dataLinkIDs;
t = nLinks + 1;
for i = 1:nLinks %create the reversed links in any case
    Link = dataLinkIDs(i,:); %link in the original FT direction
    biLink = Link; %We create a TF copy of the link with reversed attributes
    biLink(2) = Link(3);
    biLink(3) = Link(2);
    biLink(7) = Link(8);
    biLink(8) = Link(7);
    biLink(10) = mod(Link(11)+180,360);
    biLink(11) = mod(Link(10)+180,360);       
    if (Link(4)==1) %if there is a FT restriction
        %We switch the order of the TF link (->index i) and FT link (->index t) 
        %The TF link is the one that can be used and the FT link is the
        %"restricted" one
        dataFullLinksIDs(i,:)=biLink;
        dataFullLinksIDs(t,:)=Link;
    elseif (Link(4) == 2) %if there is a TF restriction
        dataFullLinksIDs(t,:)=biLink;
    else %if Link(4)=0 and there is no restriction
        dataFullLinksIDs(t,:) = biLink;
    end
        mappingLinks(i) = t;
        mappingLinks(t) = i;
        t = t + 1;       
end

% Creating incidence matrix
nStartLinks = double(NbOD);
nAbsorbingLinks = double(NbOD);
MappingDestDummy = zeros(1,nbiLinks);
A = dataFullLinksIDs(:,3); %To Nodes (end of link)
B = dataFullLinksIDs(:,2); %From Nodes (start of link)
Incidence = sparse(nbiLinks+ nStartLinks,nbiLinks + nStartLinks + nAbsorbingLinks);
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

%% Write OD demand in terms of link IDs
[row, col]=find(ODdemand>0);
ODmat=sparse(zeros(size(Incidence)));
for k=1:length(row)
    i=row(k);%origin Centroid ID
    j=col(k);%destination Centroid ID
    ODmat(nbiLinks+i,nbiLinks+nStartLinks+j)=ODdemand(i,j);
end

%% Read link couts
filename='./Input/link_counts.txt';
fid=fopen(filename);
[Nout]=textscan(fid,'%d%d%d%d%d%d%d','CommentStyle','//');
% columns: summer, spring, winter, fall, average
fclose(fid);
counts=[Nout{:}];

