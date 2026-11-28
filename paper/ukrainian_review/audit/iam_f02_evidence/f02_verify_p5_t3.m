% IAM-UA8-F02 follow-up: the preserved P5 candidate stores T3 arrays (p 26365x2, t 51512x3); compare T3 with T3 (attempt 1 compared T6).
out='<scratch>\f02';
pres='<investigator_codex_work>\2026-10-04\files-pasted-by-the-user-you\work\incremental-profile\phaseA_baseline\fresh5\step_005_candidate.mat';
P=load(pres);c=P.candidate;A=load(fullfile(out,'F02_regen_P05.mat'));F=A.F;g=F.candidate;
R=struct();
R.preservedP_vs_regenCandidateP_bitwise=isequal(c.p,g.p);R.preservedT_vs_regenCandidateT_exact=isequal(double(c.t),double(g.t));
R.preservedP_vs_meshCoord3_bitwise=isequal(c.p,F.mesh.coord3);R.preservedP_vs_localCoord3_bitwise=isequal(c.p,F.localMesh.coord3);
R.preservedT_vs_connect3_exact=isequal(double(c.t),double(F.mesh.connect3));
R.maxAbs_p_minus_regen_p_m=max(abs(c.p-g.p),[],'all');
R.crackStructEqual=isequal(c.crack,g.crack);R.gatesEqual=isequal(c.gates,g.gates);R.matEqual=isequal(c.mat,g.mat);
R.exteriorDesignEqual=isequal(c.exteriorDesign,g.exteriorDesign);R.structuredDesignEqual=isequal(c.structuredDesign,g.structuredDesign);
R.supportEqual=isequal(c.primarySupportElementIDs,g.primarySupportElementIDs);R.exteriorIDsEqual=isequal(c.exteriorElementIDs,g.exteriorElementIDs);
R.pairedNodesEqual=isequal(c.pairedNodeIDs,g.pairedNodeIDs);
R.regenFieldsMissingInPreserved=setdiff(fieldnames(g),fieldnames(c));R.preservedFieldsMissingInRegen=setdiff(fieldnames(c),fieldnames(g));
R.note='synthetic fields differ by design: the preserved candidate ran the synthetic EDI replay, the illustrative regeneration did not (RunSynthetic=false)';
fid=fopen(fullfile(out,'F02_VERIFY_P5_T3.json'),'w');fwrite(fid,jsonencode(R,'PrettyPrint',true),'char');fclose(fid);
disp(R)
