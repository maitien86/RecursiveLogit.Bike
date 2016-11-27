tripattributes=zeros(length(uniqueLinkNum),648);
for i=1:648
    for j=2:length(find(Obs(i,:)))-1
        l=linkNumber(Obs(i,j));
        tripattributes(l,i)=1;
    end
end

for i=1:648
    Names{i}=['Trip',num2str(i)];
end

nonduplicate=find(tripattributes(12386,:));
nonduplicate2=find(tripattributes(13298,:));
nonduplicate3=find(tripattributes(7907,:));
nonduplicate4=find(tripattributes(10753,:));
C=union(nonduplicate, nonduplicate2);
D=union(nonduplicate3, nonduplicate4);
E=union(C, D);
All=[1:648];
Duplicates=setdiff(All,E);