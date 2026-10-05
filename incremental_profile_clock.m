function out=incremental_profile_clock(action,varargin)
%INCREMENTAL_PROFILE_CLOCK Optional phase measurements; inactive by default.
% Explicit reset(true) enables a profiling session. Only durations and
% MATLAB memory samples are retained, never numerical or mesh state.
% Nested run/step/qualification/physical scopes have independent clocks.
persistent enabled records scopes
if isempty(enabled),enabled=false;records=empty_records();scopes=struct();end
out=[];
switch action
    case 'reset'
        validateattributes(varargin{1},{'logical'},{'scalar'});
        enabled=varargin{1};records=empty_records();scopes=struct();return
    case 'enabled',out=enabled;return
    case 'snapshot'
        out=struct2table(records);return
end
if ~enabled,return,end
scope=varargin{1};
switch action
    case 'begin'
        assert(~isfield(scopes,scope),'pathprofile:NestedScope','Scope already active.');
        scopes.(scope)=struct('segment',varargin{2},'phase','setup', ...
            'total',tic,'start',tic);
    case 'phase'
        assert(isfield(scopes,scope),'pathprofile:Scope','No active scope.');
        flush(scope,'phase');
        scopes.(scope).phase=varargin{2};scopes.(scope).start=tic;
    case 'end'
        if ~isfield(scopes,scope),return,end
        flush(scope,'phase');flush(scope,'total');scopes=rmfield(scopes,scope);
    otherwise,error('pathprofile:Action','Unknown clock action %s.',action);
end
    function flush(name,kind)
        c=scopes.(name);
        if strcmp(kind,'total'),elapsed=toc(c.total);phase='TOTAL';
        else,elapsed=toc(c.start);phase=c.phase;end
        bytes=NaN;
        if ispc
            view=memory;bytes=view.MemUsedMATLAB;
        end
        records(end+1)=struct('scope',name,'segment',c.segment,'phase',phase, ...
            'kind',kind,'seconds',elapsed,'matlabMemoryBytes',bytes);
    end
end
function r=empty_records()
r=struct('scope',{},'segment',{},'phase',{},'kind',{},'seconds',{},'matlabMemoryBytes',{});
end
