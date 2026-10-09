function assert_crack_checkpoint_physics(s,mat,C,errorID)
% Identity checks for an existing field; no change to numerical acceptance.
ok=isfield(s,'mat')&&isstruct(s.mat)&&isfield(s,'C')&&isstruct(s.C);
for name={'E','nu','ps','D','Dmat'}
    field=name{1};
    if isfield(mat,field)
        ok=ok&&isfield(s.mat,field)&&isequaln(s.mat.(field),mat.(field));
    end
end
if ok
    ok=isequaln(crack_physics_signature(s.C),crack_physics_signature(C));
end
assert(ok,errorID,'Existing checkpoint uses different frozen physical data or material.');
end
