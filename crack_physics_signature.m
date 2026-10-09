function P=crack_physics_signature(C)
% Cache identity only. Stage-II increment/mesh controls are checked separately.
% C.a0 is the historical reserved increment, not hole-only physical data.
required={'A','B','E','nu','ps','load','bc'};
P=struct();
for j=1:numel(required)
    name=required{j};assert(isfield(C,name),'crackcache:PhysicalField', ...
        'Frozen configuration is missing %s.',name);
    P.(name)=C.(name);
end
for name={'holes','hole','sig_c'}
    if isfield(C,name{1}),P.(name{1})=C.(name{1});end
end
end
