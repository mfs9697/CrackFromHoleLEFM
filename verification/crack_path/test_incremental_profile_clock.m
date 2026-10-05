function test_incremental_profile_clock()
% Instrumentation must be disabled by default and support independent
% nested scopes without including a child twice in an exclusive phase.
assert(~incremental_profile_clock('enabled'));
incremental_profile_clock('begin','disabled',0);
incremental_profile_clock('phase','disabled','unused');
incremental_profile_clock('end','disabled');
assert(isempty(incremental_profile_clock('snapshot')));
cleanup=onCleanup(@()incremental_profile_clock('reset',false)); %#ok<NASGU>
incremental_profile_clock('reset',true);
incremental_profile_clock('begin','parent',2);
incremental_profile_clock('phase','parent','outer');
incremental_profile_clock('begin','child',2);
incremental_profile_clock('phase','child','inner');
incremental_profile_clock('end','child');incremental_profile_clock('end','parent');
t=incremental_profile_clock('snapshot');
assert(all(t.seconds>=0));
assert(nnz(strcmp(t.kind,'total'))==2);
assert(all(t.segment==2));
assert(sum(t.seconds(strcmp(t.scope,'parent')&strcmp(t.kind,'phase')))<= ...
    t.seconds(strcmp(t.scope,'parent')&strcmp(t.kind,'total')));
fprintf('PASS: optional profiling clock, inactive mode and nested scopes.\n');
end
