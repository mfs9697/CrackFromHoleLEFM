function M=isolated_tip_publication_metrics(D,Ref,Vref)
% Scalar comparisons of accepted saved quantities; no extraction/refitting.
T=struct2table(D.independent);F=struct2table(D.fixed);V=[D.mouth_m(:).';T.tip_x_m,T.tip_y_m];
assert(isequal(T.segment,Ref.segment)&&size(Vref,1)>=24,'isolatedpub:ReferenceMismatch');
dx=1e6*(T.tip_x_m-Ref.tip_x_m);dy=1e6*(T.tip_y_m-Ref.tip_y_m);dr=hypot(dx,dy);
dt=1e3*(T.theta_deg-Ref.theta_deg);dq=T.KII_over_KI-Ref.KII_over_KI;
[maxDy,dyK]=max(abs(dy));[maxDr,drK]=max(dr);[maxDt,dtK]=max(abs(dt));[maxDq,dqK]=max(abs(dq));
[maxTurn,turnK]=max(abs(T.delta_theta_next_deg-Ref.delta_theta_next_deg));
[refCross,refPoint]=crossing(Ref,Vref,D.increment_m);
[isoCross,isoPoint]=crossing(T,V,D.increment_m);
M=struct('schemaVersion',1,'datasetId',D.datasetId,'authoritativeStudy',D.authoritativeStudy, ...
    'maxAbsDy_um',maxDy,'maxDySegment',dyK,'maxTipSeparation_um',maxDr,'maxTipSeparationSegment',drK, ...
    'maxAbsDtheta_mdeg',maxDt,'maxDthetaSegment',dtK,'maxAbsDq',maxDq,'maxDqSegment',dqK, ...
    'maxAbsMTSturnDifference_deg',maxTurn,'maxMTSturnDifferenceSegment',turnK, ...
    'referenceLocalSymmetry_mm',refCross,'isolatedLocalSymmetry_mm',isoCross, ...
    'localSymmetryDifference_um',1e3*(isoCross-refCross),'crossingIsInterpolated',true, ...
    'referenceCrossingPoint_mm',refPoint,'isolatedCrossingPoint_mm',isoPoint, ...
    'P23Displacement_um',[dx(23),dy(23)],'fixed',struct([]));
for scale=[2,1,.5]
    if scale==1
        rows=Ref(ismember(Ref.segment,D.fixedSegments),:);
    else
        rows=F(F.core_scale==scale,:);
    end
    j21=find(rows.segment==21);j22=find(rows.segment==22);
    zero=rows.crack_length_mm(j21)-(rows.crack_length_mm(j22)-rows.crack_length_mm(j21))* ...
        rows.KII_over_KI(j21)/(rows.KII_over_KI(j22)-rows.KII_over_KI(j21));
    j17=find(rows.segment==17);ref17=find(Ref.segment==17);
    x=struct('coreScale',scale,'q21',rows.KII_over_KI(j21),'turn21_deg',rows.delta_theta_next_deg(j21), ...
        'q22',rows.KII_over_KI(j22),'turn22_deg',rows.delta_theta_next_deg(j22),'crossing_mm',zero, ...
        'P17KIChange_percent',100*(rows.KI_unit(j17)/Ref.KI_unit(ref17)-1), ...
        'P17qDifference',rows.KII_over_KI(j17)-Ref.KII_over_KI(ref17), ...
        'P17turnDifference_deg',rows.delta_theta_next_deg(j17)-Ref.delta_theta_next_deg(ref17));
    if isempty(M.fixed),M.fixed=x;else,M.fixed(end+1)=x;end
end
M.fixedCrossingSuccessiveDifferences_um=1e3*diff([M.fixed.crossing_mm]);
M.fixedCrossingApparentOrder=log2(abs(M.fixedCrossingSuccessiveDifferences_um(1)/M.fixedCrossingSuccessiveDifferences_um(2)));
end

function [a,point]=crossing(T,V,increment)
q=T.KII_over_KI;j=find(q(1:end-1).*q(2:end)<=0,1);
assert(~isempty(j),'isolatedpub:NoCrossing');alpha=-q(j)/(q(j+1)-q(j));
a=1e3*increment*(T.segment(j)+alpha*(T.segment(j+1)-T.segment(j)));
point=1e3*((1-alpha)*V(T.segment(j)+1,:)+alpha*V(T.segment(j+1)+1,:));
end
