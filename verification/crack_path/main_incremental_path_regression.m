function Path = main_incremental_path_regression(varargin)
%MAIN_INCREMENTAL_PATH_REGRESSION
% First local regression for the general crack-path driver.
%
% This run represents three finite crack segments:
%   segment 1: theta_1=0 prescribed;
%   segment 2: theta_2 from the accepted P1 SIFs;
%   segment 3: theta_3 from the generic physical solve at P2.
%
% The generic step at P2 must reproduce the accepted Stage III-D values
% before the driver is allowed to continue to the new physical solve at P3.
%
% DEFAULT IS SAFE: AllowSolve=false.
%
% Usage:
%   Path = main_incremental_path_regression( ...
%       'FrozenState',R0,'AllowSolve',true);

    ip=inputParser;
    addParameter(ip,'FrozenState',[],@(x)isempty(x)||isstruct(x));
    addParameter(ip,'AllowSolve',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PlotEachStep',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
    parse(ip,varargin{:});
    opt=ip.Results;

    Path=run_incremental_crack_path( ...
        'FrozenState',opt.FrozenState, ...
        'MaxSegments',3, ...
        'AllowPhysicalSolves',opt.AllowSolve, ...
        'RunSynthetic',true, ...
        'ReuseCandidates',true, ...
        'RegressionGates',true, ...
        'PlotEachStep',opt.PlotEachStep, ...
        'OutputDir',opt.OutputDir);

    assert(isfield(Path.regression,'step2_pass')&&Path.regression.step2_pass, ...
        'pathreg:Stage3DRegression','Accepted Stage III-D regression failed.');
    assert(Path.nSegments==3, ...
        'pathreg:ThreeSegments','Regression must represent exactly three segments.');
    assert(height(Path.stepTable)==3, ...
        'pathreg:StepTable','Expected accepted P1 seed plus physical P2/P3 rows.');
    assert(all(Path.stepTable.pass), ...
        'pathreg:PhysicalPass','At least one physical row failed.');

    fprintf('\nGENERAL INCREMENTAL PATH REGRESSION PASS.\n');
    fprintf('  Accepted Stage III-D was reproduced at P2.\n');
    fprintf('  A genuine third segment was appended and solved at P3.\n');
    fprintf('  theta_3 = %+.12g deg\n',Path.thetaDeg(3));
    fprintf('  predicted theta_4 = %+.12g deg\n', ...
        Path.stepTable.theta_next_deg(end));
end
