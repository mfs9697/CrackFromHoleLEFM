function family=test_step62_structured_family()
%TEST_STEP62_STRUCTURED_FAMILY Geometry-only family invariants at two levels.
% No checkpoint displacement, equilibrium solve or extraction is involved.
rows=zeros(2,7);
for level=0:1
    [P,T,map,nU,rings,axisIDs,d]=build_step62_structured_patch(.006,5.4024650785e-5,level);
    [Q,S]=build_step62_structured_patch(.006,5.4024650785e-5,level);
    assert(isequal(P,Q)&&isequal(T,S),'step62:Reproducibility','Non-deterministic patch.');
    assert(nnz(any(T==1,2))==6&&numel(setdiff(unique(T(any(T==1,2),:)),1))==7);
    assert(abs(rings(1)/(5.4024650785e-5*2^-level)-1)<1e-13);
    width=diff(rings);assert(all(diff(width)>0)&&max(width(2:end)./width(1:end-1))<1.03);
    assert(all(diff(d.angularIntervalsUpper)>=0)&&max(diff(d.angularIntervalsUpper))<=3);
    assert(all(mod(d.angularIntervalsUpper,3)==0));
    neg=axisIDs(P(axisIDs,1)<0);pos=axisIDs(P(axisIDs,1)>=0);
    assert(all(map(neg)~=neg)&&all(map(pos)==pos));
    half=size(T,1)/2;U=T(1:half,:);L=T(half+1:end,[1 3 2]);
    assert(isequal(map(U),L));
    a=P(T(:,2),:)-P(T(:,1),:);b=P(T(:,3),:)-P(T(:,1),:);
    area=.5*(a(:,1).*b(:,2)-a(:,2).*b(:,1));assert(all(area>0));
    n=d.angularIntervalsUpper(end);exactArea=n*.006^2*sin(pi/n);
    assert(abs(sum(area)/exactArea-1)<1e-12);
    [P6,~]=T3toT6_fast(P,T);
    E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
    [~,~,g]=unique(E,'rows');assert(all(accumarray(g,1)<=2));
    rows(level+1,:)=[level,d.scale,size(P,1),size(T,1),size(P6,1), ...
        numel(rings),max(width(2:end)./width(1:end-1))];
end
assert(rows(2,4)>3.8*rows(1,4)&&rows(2,4)<4.2*rows(1,4));
family=array2table(rows,'VariableNames',{'level','scale','T3Nodes','T3Triangles', ...
    'T6Nodes','rings','maximumRadialWidthRatio'});
disp(family);fprintf('PASS: deterministic family, seven-edge fan, coverage and grading.\n');
end
