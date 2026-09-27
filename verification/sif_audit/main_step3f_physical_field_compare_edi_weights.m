function Out=main_step3f_physical_field_compare_edi_weights()
% Step 3F: same physical FEM field, compare analytic and FE-nodal EDI q.

here=fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));

fprintf('\n============================================================\n');
fprintf('SIF AUDIT STEP 3F: PHYSICAL-FIELD EDI WEIGHT COMPARISON\n');
fprintf('============================================================\n');

C=cfg_crack_path_two_leg_control(2.0,'plotGeom',false,'plotMesh',false);
Sol=solve_crack_path_polyline_field(C,C.sigma0);

names={"analytic_radial","fe_nodal"};
R=cell(2,1);
rows=nan(2,12);

for j=1:2
    R{j}=sweep_same_field_radii(C,C.sigma0, ...
        'radiusFractions',[0.2 0.3 0.4 0.5 0.6], ...
        'innerFactor',0.1, ...
        'innerFactorsSweep',[0.05 0.1 0.2 0.3], ...
        'fixedOuterFraction',0.5, ...
        'nthet',100, ...
        'EDIWeightFunction',names{j}, ...
        'SolvedField',Sol, ...
        'Verbose',false);

    To=R{j}.outerSweep;
    Ti=R{j}.innerSweep;

    oKIr=max(To.KI_EDI)-min(To.KI_EDI);
    oKIIr=max(To.KII_EDI)-min(To.KII_EDI);
    iKIr=max(Ti.KI_EDI)-min(Ti.KI_EDI);
    iKIIr=max(Ti.KII_EDI)-min(Ti.KII_EDI);

    [~,i05]=min(abs(To.r_over_lastLeg-0.5));

    rows(j,:)=[j, ...
        oKIr,oKIIr,iKIr,iKIIr, ...
        To.KI_old(i05),To.KI_EDI(i05), ...
        To.KII_old(i05),To.KII_EDI(i05), ...
        To.abs_dKI_over_abs_KI_EDI(i05), ...
        To.abs_dKII_over_abs_KII_EDI(i05), ...
        To.vector_difference_rel(i05)];
end

T=array2table(rows,'VariableNames',{ ...
    'weightID','outer_KI_range','outer_KII_range', ...
    'inner_KI_range','inner_KII_range', ...
    'KI_old_at_r05','KI_EDI_at_r05', ...
    'KII_old_at_r05','KII_EDI_at_r05', ...
    'rel_dKI_at_r05','rel_dKII_at_r05','vector_diff_at_r05'});

W=strings(height(T),1);
for i=1:height(T), W(i)=names{T.weightID(i)}; end
T.weightFunction=W;
T=movevars(T,'weightFunction','After','weightID');

fprintf('\nSUMMARY\n');
disp(T);

for j=1:2
    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s: matched outer-radius sweep\n',names{j});
    fprintf('------------------------------------------------------------\n');
    disp(R{j}.outerSweep(:,{ ...
        'r_over_lastLeg','KI_old','KI_EDI','KII_old','KII_EDI', ...
        'abs_dKI_over_abs_KI_EDI','abs_dKII_over_abs_KII_EDI', ...
        'vector_difference_rel'}));

    fprintf('\n%s: EDI inner-radius sweep\n',names{j});
    disp(R{j}.innerSweep);
end

Out=struct('summary',T,'runs',{R},'solution',Sol,'caseConfig',C);

fprintf('\nSTEP 3F completed.\n');
fprintf(['If fe_nodal strongly reduces EDI annulus sensitivity on the same ', ...
    'physical FEM field, promote it to the canonical EDI weight before ', ...
    'starting the deliberate mesh-asymmetry experiment.\n']);
end
