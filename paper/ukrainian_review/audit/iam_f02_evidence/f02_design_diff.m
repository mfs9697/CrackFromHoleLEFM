% IAM-UA8-F02: which design subfields differ between the preserved and the regenerated P5 candidate (read-only).
out='<scratch>\f02';
pres='<investigator_codex_work>\2026-10-04\files-pasted-by-the-user-you\work\incremental-profile\phaseA_baseline\fresh5\step_005_candidate.mat';
P=load(pres);c=P.candidate;A=load(fullfile(out,'F02_regen_P05.mat'));g=A.F.candidate;
D=struct('exteriorDesign',{diffs(c.exteriorDesign,g.exteriorDesign,'exteriorDesign')},'structuredDesign',{diffs(c.structuredDesign,g.structuredDesign,'structuredDesign')});
fid=fopen(fullfile(out,'F02_VERIFY_P5_DESIGN_DIFF.json'),'w');fwrite(fid,jsonencode(D,'PrettyPrint',true),'char');fclose(fid);
disp(D.exteriorDesign);disp(D.structuredDesign);
function L=diffs(a,b,p)
L={};
if isstruct(a)&&isstruct(b)&&isscalar(a)&&isscalar(b)
 fa=fieldnames(a);fb=fieldnames(b);
 for f=reshape(setdiff(fa,fb),1,[]),L{end+1}=[p '.' f{1} ' (only preserved)'];end %#ok<AGROW>
 for f=reshape(setdiff(fb,fa),1,[]),L{end+1}=[p '.' f{1} ' (only regenerated)'];end %#ok<AGROW>
 for f=reshape(intersect(fa,fb),1,[]),L=[L,diffs(a.(f{1}),b.(f{1}),[p '.' f{1}])];end %#ok<AGROW>
elseif ~isequaln(a,b)
 if isnumeric(a)&&isnumeric(b)&&isequal(size(a),size(b)),L{end+1}=sprintf('%s maxAbsDiff=%g',p,max(abs(double(a(:))-double(b(:)))));
 elseif islogical(a)||islogical(b),L{end+1}=sprintf('%s preserved=%s regenerated=%s',p,mat2str(a),mat2str(b));
 else,L{end+1}=[p ' differs'];end
end
end
