function [D,T,F]=load_isolated_tip_publication_data(file)
% Strict portable isolated-core reader. No historical archive/CSV fallback.
assert(exist(file,'file')==2,'isolatedpub:MissingDataset','Missing isolated publication dataset: %s',file);
folder=fileparts(file);manifestFile=fullfile(folder,'isolated_tip_resolution_manifest.json');
assert(exist(manifestFile,'file')==2,'isolatedpub:MissingManifest','Missing isolated provenance manifest.');
D=jsondecode(fileread(file));M=jsondecode(fileread(manifestFile));
assert(D.schemaVersion==1&&M.schemaVersion==1&&strcmp(D.datasetId,'isolated_core_tip_resolution')&& ...
    strcmp(M.datasetId,D.datasetId)&&strcmp(D.authoritativeStudy,M.authoritativeStudy)&& ...
    D.exteriorScale==1&&D.requestedExteriorLawFixed&&~D.identicalExteriorConnectivityRequired, ...
    'isolatedpub:FamilyMismatch','Dataset is not the isolated core / reference requested exterior law.');
for item=M.portableFiles(:).'
    [~,name,extension]=fileparts(item.path);input=fullfile(folder,[name extension]);
    assert(exist(input,'file')==2&&strcmp(publication_sha256(input),item.sha256), ...
        'isolatedpub:HashMismatch','Portable data bytes differ from manifest: %s',input);
end
assert(any(arrayfun(@(x)strcmp(x.sha256,publication_sha256(file)),M.portableFiles)), ...
    'isolatedpub:HashMismatch','Dataset itself is not covered by manifest.');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
for item=M.referenceInputs(:).'
    input=fullfile(root,strrep(item.path,'/',filesep));
    assert(exist(input,'file')==2&&strcmp(publication_sha256(input),item.sha256), ...
        'isolatedpub:ReferenceMismatch','Reference evidence differs from the source manifest.');
end
assert(D.independentPhysicalStates==23&&D.referencePhysicalStates==23&& ...
    isequal(D.fixedSegments(:),[17;21;22;23]),'isolatedpub:Incomplete','Incomplete publication study.');
T=struct2table(D.independent);F=struct2table(D.fixed);
assert(height(T)==23&&isequal(T.segment,(1:23)')&&all(T.core_scale==2)&&all(T.exterior_scale==1), ...
    'isolatedpub:FamilyMismatch','Independent history is not complete core=2/exterior=1.');
assert(height(F)==8&&all(F.exterior_scale==1),'isolatedpub:Incomplete','Eight fixed isolated results required.');
for scale=[2,.5]
    assert(isequal(F.segment(F.core_scale==scale),[17;21;22;23]), ...
        'isolatedpub:FixedStates','Wrong fixed-reference geometry states.');
end
for rows={D.independent,D.fixed}
    for r=rows{1}(:).'
        assert(r.pass&&r.qualificationPassed&&r.syntheticPassed&& ...
            all(structfun(@logical,r.physicalGates))&&r.pcg_flag==0&&isfinite(r.true_rel_residual), ...
            'isolatedpub:Acceptance','An unaccepted physical record cannot be plotted.');
        assert(abs(r.KII_over_KI-r.KII_unit/r.KI_unit)<=1e-14&&r.KI_unit>0, ...
            'isolatedpub:RecordMismatch','Inconsistent saved SIF ratio.');
        if r.core_scale==2,expected=[19;28;23;18];else,expected=[74;108;86;67];end
        assert(isequal(r.nativePoints(:),expected),'isolatedpub:RecordMismatch','Wrong native sampling fingerprint.');
    end
end
V=[D.mouth_m(:).';T.tip_x_m,T.tip_y_m];
assert(all(abs(vecnorm(diff(V),2,2)-D.increment_m)<=2e-12)&& ...
    all(abs(T.crack_length_mm-1e3*D.increment_m*T.segment)<=1e-8), ...
    'isolatedpub:RecordMismatch','Trajectory lengths do not match physical states.');
end
