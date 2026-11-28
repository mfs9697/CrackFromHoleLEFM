function R = iam_f02_regenerate_mesh(genRoot, outDir, preservedP5)
%IAM_F02_REGENERATE_MESH  Deterministic ILLUSTRATIVE regeneration of the accepted M0 mesh P21 for the proposed Figure 2.
%   R = iam_f02_regenerate_mesh(genRoot, outDir)                 regenerate P21, verify the production fingerprint, export
%   R = iam_f02_regenerate_mesh(genRoot, outDir, preservedP5)    also regenerate P5 and compare it bitwise with a preserved
%                                                                historical P5 candidate (step_005_candidate.mat)
% genRoot  : clean checkout of the pinned historical generator cfdf0f110010d688a4e3c48f6d88a00fd17dc698 (its pinned sources equal
%            paper/data/source_manifest.csv in CRLF form); refused otherwise.
% outDir   : fresh folder for the raw little-endian float64 exports (RAW_INDEX.json) and IAM_F02_REGENERATION.json.
% Inputs   : paper/data/evidence_exact.mat (accepted vertices P0..P24, full-precision qualification records) and
%            paper/data/accepted_stage1_source.mat (frozen Stage-I R0) of THIS repository (read-only).
% Only mesh construction and its structural gates: qualify_incremental_crack_candidate with RunSynthetic=false (no synthetic
% EDI replay), SaveCandidate=false, SaveCompact=false. No stiffness assembly, no FEM solve, no SIF/EDI evaluation.
% The production meshes themselves were never archived (final_clean_run kept compact records only); the regenerated mesh is an
% illustration whose identity with the production mesh is supported by fingerprints, not proven bitwise.
    here = fileparts(mfilename('fullpath')); repo = fileparts(fileparts(fileparts(here)));
    PIN = 'cfdf0f110010d688a4e3c48f6d88a00fd17dc698';
    [s1, h] = system(sprintf('git -C "%s" rev-parse HEAD', genRoot)); [s2, st] = system(sprintf('git -C "%s" status --porcelain', genRoot));
    assert(s1 == 0 && strcmp(strtrim(h), PIN) && s2 == 0 && isempty(strtrim(st)), 'iamf02:Generator', 'genRoot must be a clean checkout of %s.', PIN);
    assert(~isfolder(outDir) || isempty(dir(fullfile(outDir, '*.f64'))), 'iamf02:Out', 'outDir must be fresh.');
    if ~isfolder(outDir), mkdir(outDir); end
    oldPath = path; cleanup = onCleanup(@() path(oldPath)); restoredefaultpath; addpath(genpath(genRoot));
    assert(startsWith(which('qualify_incremental_crack_candidate'), genRoot), 'iamf02:Path', 'Generator not first on the path.');
    Z = load(fullfile(repo, 'paper', 'data', 'accepted_stage1_source.mat'), 'R0'); R0 = Z.R0;
    Ev = load(fullfile(repo, 'paper', 'data', 'evidence_exact.mat'), 'E'); E = Ev.E; V = E.vertices_m;
    Q = E.qualification; if isstruct(Q), Q = struct2table(Q); end
    cols = {'T3_nodes','T3_elements','T6_nodes','core_T3_elements','exterior_T3_elements','EDI_elements','physical_clearance_m','exterior_crack_m', ...
        'area_rel_error','min_angle_deg','max_neighbor_ratio','min_T6_detJ_m2','core_coord_error_m','pair_error_T3_m','pair_error_T6_m'};
    R = struct('generator', PIN, 'noFEM', true, 'runSynthetic', false);
    states = 21; if nargin > 2 && ~isempty(preservedP5), states = [5 21]; end
    for k = states
        F = qualify_incremental_crack_candidate(V(1:k+1, :), 'FrozenState', R0, 'RunSynthetic', false, 'SaveCandidate', false, ...
            'SaveCompact', false, 'Plot', false, 'ExteriorVerbose', false);
        S = F.summary; row = Q(Q.segment == k, :); fp = struct();
        for j = 1:numel(cols), fp.(cols{j}) = struct('regenerated', S.(cols{j}), 'production', row.(cols{j}), 'equal', isequal(S.(cols{j}), row.(cols{j}))); end
        e = struct('structuralGatesPass', all(structfun(@logical, F.gates)), 'fingerprint', fp, ...
            'fingerprintAllEqual', all(cellfun(@(f) fp.(f).equal, cols)));
        if k == 5
            P = load(preservedP5); c = P.candidate; g = F.candidate;
            e.preservedPathEqualsAccepted = isequal(c.path, V(1:6, :));
            e.T3NodesBitwise = isequal(c.p, F.mesh.coord3); e.T3ConnectivityExact = isequal(double(c.t), double(F.mesh.connect3));
            e.crackEqual = isequal(c.crack, g.crack); e.gatesEqual = isequal(c.gates, g.gates); e.pairedNodesEqual = isequal(c.pairedNodeIDs, g.pairedNodeIDs);
            e.mirrorMapEqual = isequal(c.pairedMirrorLocal, g.pairedMirrorLocal); e.EDIsupportEqual = isequal(c.primarySupportElementIDs, g.primarySupportElementIDs);
        else
            c = F.candidate; IX = struct('file', {}, 'shape', {});
            IX = ex(IX, outDir, 'X3_mm', 1e3*F.mesh.coord3); IX = ex(IX, outDir, 'T3', double(F.mesh.connect3)); IX = ex(IX, outDir, 'XL3_mm', 1e3*F.localMesh.coord3);
            IX = ex(IX, outDir, 'upperNodes', double(F.crackLocal.upperNodes(:))); IX = ex(IX, outDir, 'lowerNodes', double(F.crackLocal.lowerNodes(:)));
            IX = ex(IX, outDir, 'tipNode', double(F.crackLocal.tipNode)); IX = ex(IX, outDir, 'pathGlobal_mm', 1e3*c.path);
            IX = ex(IX, outDir, 'pairedElementIDs', double(c.pairedElementIDs(:))); IX = ex(IX, outDir, 'exteriorElementIDs', double(c.exteriorElementIDs(:)));
            IX = ex(IX, outDir, 'EDIsupportIDs', double(c.primarySupportElementIDs(:))); IX = ex(IX, outDir, 'pairedNodeIDs', double(c.pairedNodeIDs(:)));
            IX = ex(IX, outDir, 'pairedMirrorLocal', double(c.pairedMirrorLocal));
            fid = fopen(fullfile(outDir, 'RAW_INDEX.json'), 'w'); fwrite(fid, jsonencode(IX), 'char'); fclose(fid);
        end
        R.(sprintf('P%d', k)) = e;
        fprintf('IAM_F02 P%d gates=%d fingerprint=%d\n', k, e.structuralGatesPass, e.fingerprintAllEqual);
    end
    fid = fopen(fullfile(outDir, 'IAM_F02_REGENERATION.json'), 'w'); fwrite(fid, jsonencode(R, 'PrettyPrint', true), 'char'); fclose(fid);
end

function IX = ex(IX, d, name, x)
    f = fullfile(d, [name '.f64']); fid = fopen(f, 'w', 'l'); fwrite(fid, double(x(:)), 'float64'); fclose(fid);
    IX(end+1) = struct('file', [name '.f64'], 'shape', size(x));
end
