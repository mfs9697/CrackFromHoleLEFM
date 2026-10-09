function Report=test_isolated_tip_publication(varargin)
% Publication-only guards, including a checkout containing no local archives.
ip=inputParser;addParameter(ip,'WorkDir',tempname,@(x)ischar(x)||isstring(x));parse(ip,varargin{:});
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));work=char(ip.Results.WorkDir);mkdir(work);
folder=fullfile(root,'paper','data');file=fullfile(folder,'isolated_tip_resolution.json');
[D,T,F]=load_isolated_tip_publication_data(file);assert(height(T)==23&&height(F)==8);count=1;
z=load(fullfile(folder,'evidence_exact.mat'),'E','T');M=isolated_tip_publication_metrics(D,z.T,z.E.vertices_m);
assert(abs(M.maxTipSeparation_um-.127996714956727)<1e-12&& ...
    abs(M.isolatedLocalSymmetry_mm-84.14627637907)<1e-10);count=count+1;
expect(@()plot_tip2h0_vs_reference_publication('IsolatedDatasetFile',fullfile(work,'absent.json'),'Export',false), ...
    'isolatedpub:MissingDataset');count=count+1;
expect(@()plot_tip2h0_vs_reference_publication('Tip2h0StateFile', ...
    fullfile(root,'verification','crack_path','tip_2h0_independent_run','trajectory','path_run_state.mat'),'Export',false), ...
    'tip2h0fig:HistoricalInputRejected');count=count+1;
badDir=fullfile(work,'mutants');mkdir(badDir);
for name={'isolated_tip_resolution.json','isolated_tip_resolution_manifest.json','isolated_tip_fixed.csv','isolated_tip_states.csv'}
    copyfile(fullfile(folder,name{1}),fullfile(badDir,name{1}));
end
badFile=fullfile(badDir,'isolated_tip_resolution.json');
bad=D;bad.independent(1).exterior_scale=2;write_json(badFile,bad);
expect(@()load_isolated_tip_publication_data(badFile),'isolatedpub:HashMismatch');count=count+1;
rehash();expect(@()load_isolated_tip_publication_data(badFile),'isolatedpub:FamilyMismatch');count=count+1;
bad=D;bad.independent=bad.independent(1:22);write_json(badFile,bad);rehash();
expect(@()load_isolated_tip_publication_data(badFile),'isolatedpub:FamilyMismatch');count=count+1;
bad=D;keys=fieldnames(bad.independent(1).physicalGates);bad.independent(1).physicalGates.(keys{1})=false;
write_json(badFile,bad);rehash();expect(@()load_isolated_tip_publication_data(badFile),'isolatedpub:Acceptance');count=count+1;
bad=D;bad.datasetId='historical_coupled_tip_family';write_json(badFile,bad);rehash();
expect(@()load_isolated_tip_publication_data(badFile),'isolatedpub:FamilyMismatch');count=count+1;
Report=struct('pass',true,'checks',count,'physicalSolves',0);
save(fullfile(work,'publication_test_report.mat'),'Report','-v7');
fprintf('ISOLATED PUBLICATION TESTS PASS: %d checks; no physical calculations.\n',count);

    function rehash()
        mf=fullfile(badDir,'isolated_tip_resolution_manifest.json');m=jsondecode(fileread(mf));
        for j=1:numel(m.portableFiles)
            [~,name,ext]=fileparts(m.portableFiles(j).path);
            if strcmp([name ext],'isolated_tip_resolution.json'),m.portableFiles(j).sha256=publication_sha256(badFile);end
        end
        write_json(mf,m);
    end
end
function expect(fn,id)
try,fn();catch ME
    assert(strcmp(ME.identifier,id),'isolatedpub:TestError','Expected %s, got %s: %s',id,ME.identifier,ME.message);return
end
error('isolatedpub:MissingExpectedError','Expected %s.',id);
end
function write_json(file,s)
fid=fopen(file,'wb');cleanup=onCleanup(@()fclose(fid));fwrite(fid,unicode2native(jsonencode(s),'UTF-8'),'uint8');
end
