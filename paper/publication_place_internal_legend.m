function info=publication_place_internal_legend(ax,lg)
% Graphics-only collision check; never alters XData/YData or scientific gates.
ax.Units='centimeters';A=ax.Position;lg.Units='centimeters';drawnow;
sz=lg.Position(3:4);xl=ax.XLim;yl=ax.YLim;
if strcmp(ax.DataAspectRatioMode,'manual')
    aspect=(diff(xl)/ax.DataAspectRatio(1))/(diff(yl)/ax.DataAspectRatio(2));
    if A(3)/A(4)>aspect,w=A(4)*aspect;A(1)=A(1)+(A(3)-w)/2;A(3)=w;
    else,h=A(3)/aspect;A(2)=A(2)+(A(4)-h)/2;A(4)=h;end
end
xy=@(x,y)[A(1)+(x(:)-xl(1))/diff(xl)*A(3),A(2)+(y(:)-yl(1))/diff(yl)*A(4)];
paths={};reserved={};
for obj=findall(ax,'Type','line').'
    p=xy(obj.XData,obj.YData);p=p(all(isfinite(p),2),:);
    if isempty(p),continue;end
    paths{end+1}=p; %#ok<AGROW>
    if size(p,1)>30&&norm(p(1,:)-p(end,:))<.01
        reserved{end+1}=[min(p(:,1)),min(p(:,2)),max(p(:,1)),max(p(:,2))]; %#ok<AGROW>
    end
end
for obj=findall(ax,'Type','constantline').'
    if strcmp(obj.InterceptAxis,'x'),p=xy([obj.Value obj.Value],yl);
    else,p=xy(xl,[obj.Value obj.Value]);end
    paths{end+1}=p; %#ok<AGROW>
end
for obj=findall(ax,'Type','scatter').'
    p=xy(obj.XData,obj.YData);p=p(all(isfinite(p),2),:);
    for k=1:size(p,1),reserved{end+1}=[p(k,:)-.06,p(k,:)+.06];end %#ok<AGROW>
end
for obj=findall(ax,'Type','text').'
    if ismember(obj,[ax.XLabel ax.YLabel ax.Title]),continue;end
    old=obj.Units;obj.Units='data';e=obj.Extent;obj.Units=old;
    p=xy([e(1),e(1)+e(3)],[e(2),e(2)+e(4)]);reserved{end+1}=[p(1,:),p(2,:)]; %#ok<AGROW>
end
names={'northeast','northwest','southwest','southeast'};
for pad=[.15 .25 .30 .40]
    positions={[A(1)+A(3)-sz(1)-pad,A(2)+A(4)-sz(2)-pad], ...
        [A(1)+pad,A(2)+A(4)-sz(2)-pad],[A(1)+pad,A(2)+pad], ...
        [A(1)+A(3)-sz(1)-pad,A(2)+pad]};
    for j=1:4
        p=positions{j};r=[p,p+sz];
        if r(1)<A(1)||r(2)<A(2)||r(3)>A(1)+A(3)||r(4)>A(2)+A(4),continue;end
        test=r+[-.025 -.025 .025 .025];hit=false;
        for k=1:numel(reserved)
            b=reserved{k};if test(1)<b(3)&&b(1)<test(3)&&test(2)<b(4)&&b(2)<test(4),hit=true;break;end
        end
        if hit,continue;end
        for k=1:numel(paths)
            q=paths{k};
            if size(q,1)==1,hit=all(q>=test(1:2))&&all(q<=test(3:4));end
            for m=1:size(q,1)-1
                if segment_intersects(q(m,:),q(m+1,:),test),hit=true;break;end
            end
            if hit,break;end
        end
        if ~hit
            lg.Position=[p,sz];info=struct('bounds_mm',lg.Position*10, ...
                'insideAxes',true,'collisionCount',0,'location',names{j});return
        end
    end
end
error('publication:LegendCollision','No unobstructed internal legend position; inspect the review instead of hiding data.');
end
function hit=segment_intersects(p,q,r)
d=q-p;lo=0;hi=1;aa=[-d(1),d(1),-d(2),d(2)];bb=[p(1)-r(1),r(3)-p(1),p(2)-r(2),r(4)-p(2)];
for k=1:4
    if abs(aa(k))<1e-12,if bb(k)<0,hit=false;return;end
    elseif aa(k)<0,lo=max(lo,bb(k)/aa(k));else,hi=min(hi,bb(k)/aa(k));end
    if lo>hi,hit=false;return;end
end
hit=true;
end
