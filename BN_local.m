function [B, DetJ, dNdx] = BN_local(xi0, X)
%BN_LOCAL  T6 strain-displacement matrix and shape gradients.
%
% Node ordering:
%   1=L1 vertex, 2=L2 vertex, 3=L3 vertex,
%   4=edge 1-2, 5=edge 2-3, 6=edge 3-1.
%
% Inputs
%   xi0  [2x1] parent coordinates [L1; L2], L3=1-L1-L2
%   X    [6x2] nodal coordinates
%
% Outputs
%   B      [3x12] engineering-strain matrix [exx; eyy; gxy]
%   DetJ   Jacobian determinant
%   dNdx   [2x6] derivatives wrt physical x,y

    xi0 = xi0(:);
    if numel(xi0) ~= 2
        error('BN_local:BadXi', 'xi0 must contain [L1; L2].');
    end
    if ~isequal(size(X), [6,2])
        error('BN_local:BadX', 'X must be 6x2 for a T6 triangle.');
    end

    L1 = xi0(1);
    L2 = xi0(2);
    L3 = 1 - L1 - L2;

    dN_dL1 = [ ...
        4*L1 - 1, ...
        0, ...
        -(4*L3 - 1), ...
        4*L2, ...
        -4*L2, ...
        4*(L3 - L1)];

    dN_dL2 = [ ...
        0, ...
        4*L2 - 1, ...
        -(4*L3 - 1), ...
        4*L1, ...
        4*(L3 - L2), ...
        -4*L1];

    dNdxi = [dN_dL1; dN_dL2];

    J = dNdxi * X;
    DetJ = det(J);

    if abs(DetJ) <= eps(max(1, norm(J, 'fro')))
        error('BN_local:DegenerateElement', ...
            'Degenerate T6 mapping (DetJ = %.3e).', DetJ);
    end

    dNdx = J \ dNdxi;

    inx = 1:2:11;
    iny = 2:2:12;

    B = zeros(3,12);
    B(1,inx) = dNdx(1,:);
    B(2,iny) = dNdx(2,:);
    B(3,inx) = dNdx(2,:);
    B(3,iny) = dNdx(1,:);
end
