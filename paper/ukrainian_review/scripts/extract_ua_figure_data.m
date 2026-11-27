function M=extract_ua_figure_data(varargin)
%EXTRACT_UA_FIGURE_DATA  Plotting tables for the Ukrainian figure rebuild (AAM-UA7). Read-only on every record.
%  M = extract_ua_figure_data()                     writes paper/ukrainian_review/figures_src/data/*.csv + DATA_MANIFEST.json
%  M = extract_ua_figure_data('OutDir',dir)         another (empty) output folder, e.g. for negative controls
% Sources (authoritative, unchanged): paper/data/evidence_exact.mat (E,T,COD,Q; same loader as the English
% plot_existing_manuscript_figures.m), paper/data/isolated_tip_resolution.json (load_isolated_tip_publication_data),
% paper/data/increment_plot_records.mat (Increment; same fields as plot_saved_increment_records.m).
% No FE solve, mesh generation, EDI replay, COD refit, MTS recomputation, interpolation, smoothing or extrapolation.
% Transformations are explicit linear unit/scale factors only (listed in the manifest): m->mm (1e3), m->um (1e6),
% deg->mdeg (1e3), dimensionless scale factors 1e3 / 1e5 for displayed ratios, residual / 5e-10 (gate stated in the
% manuscript). Values are written with %.17g (round-trip exact for IEEE double).
ip=inputParser;addParameter(ip,'OutDir','',@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
here=fileparts(mfilename('fullpath'));ua=fileparts(here);paper=fileparts(ua);root=fileparts(paper);
out=char(ip.Results.OutDir);if isempty(out),out=fullfile(ua,'figures_src','data');end
if ~isfolder(out),mkdir(out);end
assert(isempty(dir(fullfile(out,'*.csv'))),'uafig:Overwrite','Output folder already contains CSV tables: %s',out);
src=struct();
ef=fullfile(paper,'data','evidence_exact.mat');jf=fullfile(paper,'data','isolated_tip_resolution.json');incf=fullfile(paper,'data','increment_plot_records.mat');
for f={ef,jf,incf},assert(isfile(f{1}),'uafig:Missing','Missing record %s',f{1});end
A=load(ef,'E','T','COD','Q');E=A.E;T=A.T;COD=A.COD;Q=A.Q;
assert(height(T)==23&&isequal(T.segment(:),(1:23)')&&height(COD)==176&&E.modeMixityPeakSegment==17,'uafig:Evidence','Unexpected evidence structure.');
V=E.vertices_m;assert(isequal(size(V),[25 2]),'uafig:Vertices','Expected P0..P24 (25 vertices).');
M=struct('schema','AAM-UA7.FigureData.v1','sources',struct(),'tables',struct(),'constants',struct());
M.sources.evidence=fileInfo(ef,root);M.sources.isolated=fileInfo(jf,root);M.sources.increment=fileInfo(incf,root);
% ---------------------------------------------------------------- constants (geometry and gates)
hole=E.C.hole;c=struct('hole_cx_mm',1e3*hole.center(1),'hole_cy_mm',1e3*hole.center(2),'hole_r_mm',1e3*hole.r, ...
 'right_boundary_x_mm',1e3*E.C.A,'aLS_mm',E.interpolatedZeroLength_mm,'xLS_mm',1e3*E.interpolatedZeroPoint_m(1),'yLS_mm',1e3*E.interpolatedZeroPoint_m(2), ...
 'q_peak_segment',E.modeMixityPeakSegment,'true_residual_gate',5e-10,'min_angle_gate_deg',20,'size_ratio_gate',1.8,'increment_mm',1e3*E.increment_m);
M.constants=c;
% ---------------------------------------------------------------- Fig. 3 trajectory (P0..P23 accepted, P23-P24 unsolved)
k=(0:24)';status=[repmat({'accepted'},24,1);{'unsolved'}];
t=table(k,1e3*V(:,1),1e3*V(:,2),status,'VariableNames',{'k','x_mm','y_mm','status'});
wt('fig03_vertices',t,'k = P_k index; x_mm=1e3*E.vertices_m(:,1), y_mm=1e3*E.vertices_m(:,2); P24 = qualified unsolved extension end (no accepted solution)');
% ---------------------------------------------------------------- Fig. 4 KI, KII, theta (23 accepted states)
t=table(T.segment,T.crack_length_mm,T.KI_unit,T.KII_unit,1e3*T.KII_unit,T.theta_deg,'VariableNames',{'segment','a_mm','KI','KII','KII_1e3','theta_deg'});
wt('fig04_states',t,'T.crack_length_mm, T.KI_unit, T.KII_unit (MPa*m^1/2 at unit load), KII_1e3 = 1e3*T.KII_unit, T.theta_deg');
% ---------------------------------------------------------------- Fig. 5 late mode mixity and MTS turn (a >= 60 mm)
L=T(T.crack_length_mm>=60,:);
t=table(L.segment,L.crack_length_mm,L.KII_over_KI,1e3*L.KII_over_KI,L.delta_theta_next_deg,'VariableNames',{'segment','a_mm','q','q_1e3','dtheta_next_deg'});
wt('fig05_late',t,'rows with T.crack_length_mm >= 60; q_1e3 = 1e3*T.KII_over_KI; dtheta_next_deg = T.delta_theta_next_deg');
% ---------------------------------------------------------------- Fig. 6 COD-EDI differences (a >= 60 mm, 4 windows x 2 degrees)
D=COD(COD.crack_length_mm>=60,:);W=[0.04 0.20;0.04 0.30;0.08 0.30;0.12 0.30];rows={};
for deg=1:2
 for j=1:4
  u=D.degree==deg&abs(D.lower_r_over_DeltaA-W(j,1))<1e-12&abs(D.upper_r_over_DeltaA-W(j,2))<1e-12;
  X=D.crack_length_mm(u);[X,o]=sort(X);te=D.turn_error_deg(u);te=te(o);re=D.ratio_error(u);re=re(o);sg=D.segment(u);sg=sg(o);
  for i=1:numel(X),rows(end+1,:)={deg,j,W(j,1),W(j,2),sg(i),X(i),te(i),1e3*te(i),re(i),1e5*re(i)};end %#ok<AGROW>
 end
end
t=cell2table(rows,'VariableNames',{'degree','window','lower','upper','segment','a_mm','turn_err_deg','turn_err_1e3deg','ratio_err','ratio_err_1e5'});
wt('fig06_cod',t,'COD rows with crack_length_mm >= 60, selected by degree and window (|bound-W|<1e-12), sorted by a; turn_err_1e3deg = 1e3*turn_error_deg; ratio_err_1e5 = 1e5*ratio_error');
% ---------------------------------------------------------------- Fig. 7 tip-resolution sensitivity (h0 = reference, 2h0 = isolated core)
J=jsondecode(fileread(jf));I=struct2table(J.independent);I=sortrows(I,'segment');
assert(height(I)==23&&isequal(I.segment(:),(1:23)')&&all(I.core_scale==2)&&all(I.exterior_scale==1),'uafig:Isolated','Expected 23 independent 2h0 states (core scale 2, exterior 1).');
Vm=[J.mouth_m(:).';[I.tip_x_m,I.tip_y_m]];Vr=V(1:24,:);
assert(isequal(Vm(1,:),Vr(1,:)),'uafig:Mouth','Mouth points differ.');
t=table((0:23)',1e3*Vr(:,1),1e3*Vr(:,2),1e3*Vm(:,1),1e3*Vm(:,2),'VariableNames',{'k','x_h0_mm','y_h0_mm','x_2h0_mm','y_2h0_mm'});
wt('fig07_vertices',t,'h0: 1e3*E.vertices_m(1:24,:); 2h0: 1e3*[mouth_m; independent tip_x_m, tip_y_m] (isolated_tip_resolution.json)');
da=norm(Vr(2,:)-Vr(1,:));a=1e3*da*T.segment;
t=table(T.segment,a,1e6*(I.tip_y_m-T.tip_y_m),1e3*(I.theta_deg-T.theta_deg),T.KII_over_KI,I.KII_over_KI,1e3*T.KII_over_KI,1e3*I.KII_over_KI, ...
 'VariableNames',{'segment','a_mm','dy_um','dtheta_mdeg','q_h0','q_2h0','q_h0_1e3','q_2h0_1e3'});
wt('fig07_metrics',t,'a_mm = 1e3*|V1-V0|*segment (as plot_tip2h0_vs_reference_publication.m); dy_um = 1e6*(y_2h0 - y_h0); dtheta_mdeg = 1e3*(theta_2h0 - theta_h0); q = KII/KI; q_*_1e3 = 1e3*q');
[~,ir]=max(T.KII_over_KI);[~,im]=max(I.KII_over_KI);
c2=struct('qmax_h0_segment',T.segment(ir),'qmax_h0_x_mm',1e3*T.tip_x_m(ir),'qmax_h0_y_mm',1e3*T.tip_y_m(ir), ...
 'qmax_2h0_segment',I.segment(im),'qmax_2h0_x_mm',1e3*I.tip_x_m(im),'qmax_2h0_y_mm',1e3*I.tip_y_m(im),'last_h0_x_mm',1e3*Vr(end,1),'last_h0_y_mm',1e3*Vr(end,2));
M.constants.fig07=c2;
% ---------------------------------------------------------------- Fig. 9 increment companions (saved plotting records)
z=load(incf,'Increment');Ic=z.Increment;Tm=Ic.metrics;N=Ic.native;
assert(Ic.has1mm&&Ic.maxCommonCrackLength_mm==92&&height(Tm)==23&&height(N)==161,'uafig:Increment','Incomplete increment records.');
t=table(Tm.crack_length_mm,Tm.dy_2minus4_um,Tm.dy_1minus2_um,Tm.dtheta_2minus4_mdeg,Tm.dtheta_1minus2_mdeg, ...
 'VariableNames',{'a_mm','dy_2m4_um','dy_1m2_um','dth_2m4_mdeg','dth_1m2_mdeg'});
wt('fig09_metrics',t,'Increment.metrics: crack_length_mm, dy_2minus4_um, dy_1minus2_um, dtheta_2minus4_mdeg, dtheta_1minus2_mdeg; no transformation');
t=table(N.increment_mm,N.crack_length_mm,N.turn_density_deg_per_mm,'VariableNames',{'increment_mm','a_mm','turn_density_deg_per_mm'});
wt('fig09_native',t,'Increment.native rows (record order): increment_mm, crack_length_mm, turn_density_deg_per_mm; no transformation');
% ---------------------------------------------------------------- Fig. 10 numerical quality
Qp=Q(Q.segment<=23,:);Tp=T(T.segment>=2,:);
t=table(Qp.segment,1e3*Qp.path_length_m,Qp.min_angle_deg,Qp.max_neighbor_ratio,'VariableNames',{'segment','s_mm','min_angle_deg','max_neighbor_ratio'});
wt('fig10_mesh',t,'Q rows with segment <= 23: s_mm = 1e3*Q.path_length_m, Q.min_angle_deg, Q.max_neighbor_ratio');
t=table(Tp.segment,Tp.crack_length_mm,Tp.PCG_iterations,Tp.true_rel_residual,Tp.true_rel_residual/5e-10, ...
 'VariableNames',{'segment','a_mm','pcg_iterations','true_rel_residual','residual_over_gate'});
wt('fig10_solver',t,'T rows with segment >= 2 (newly solved states): PCG_iterations, true_rel_residual, residual_over_gate = true_rel_residual/5e-10 (gate stated in the manuscript)');
fid=fopen(fullfile(out,'DATA_MANIFEST.json'),'w');fwrite(fid,jsonencode(M,'PrettyPrint',true),'char');fclose(fid);
fprintf('UA7_EXTRACT tables=%d out=%s\n',numel(fieldnames(M.tables)),out);

 function wt(name,t,how)
  f=fullfile(out,[name '.csv']);fid=fopen(f,'w');fprintf(fid,'%s\n',strjoin(t.Properties.VariableNames,','));
  for r=1:height(t)
   parts=cell(1,width(t));
   for q=1:width(t),v=t{r,q};if iscell(v),v=v{1};end
    if ischar(v)||isstring(v),parts{q}=char(v);else,parts{q}=sprintf('%.17g',v);end
   end
   fprintf(fid,'%s\n',strjoin(parts,','));
  end
  fclose(fid);M.tables.(name)=struct('file',[name '.csv'],'rows',height(t),'columns',{t.Properties.VariableNames},'transformation',how, ...
   'sha256',lower(g2sha(f)));
 end
end

function s=fileInfo(f,root)
d=dir(f);s=struct('path',strrep(strrep(f,[root filesep],''),'\','/'),'sha256',lower(g2sha(f)),'bytes',d.bytes);
end
function h=g2sha(f)
fid=fopen(f,'r');b=fread(fid,inf,'*uint8');fclose(fid);md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
h=reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]);
end
