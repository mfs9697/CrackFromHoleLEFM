function [B,Det,dNdx] = BN_local(xi0,X)
%BN_LOCAL T6 strain-displacement matrix and shape-function gradients.
%
%   [B,Det] = BN_local(xi0,X)
%   [B,Det,dNdx] = BN_local(xi0,X)
%
% Inputs
%   xi0 : 2-by-1 parent coordinates [xi; eta]. The third barycentric
%         coordinate is 1-xi-eta.
%   X   : 6-by-2 coordinates of one quadratic triangular (T6) element,
%         ordered [1 2 3 12 23 31].
%
% Outputs
%   B     : 3-by-12 engineering-strain matrix
%   Det   : determinant of the isoparametric Jacobian
%   dNdx  : 2-by-6 global shape-function gradients
%           row 1 = dN/dx, row 2 = dN/dy
%
% This shared helper is the repository-level version of the historical
% local helper embedded in StressExt.m. Several physical/SIF routines call
% BN_local as an external function, so keeping it here makes a fresh clone
% self-contained. The algebra is intentionally identical to the audited
% historical helper; the only extension is the optional third output dNdx.

    validateattributes(xi0,{'numeric'},{'vector','numel',2,'finite','real'});
    validateattributes(X,{'numeric'},{'size',[6 2],'finite','real'});

    xi0=xi0(:);
    xi=[xi0;1-sum(xi0)];

    % Parent-coordinate derivatives for the T6 shape functions.
    Nap=[ ...
        4*xi(1)-1, 0,           1-4*xi(3), ...
        4*xi(2),  -4*xi(2),     4*xi(3)-4*xi(1); ...
        0,         4*xi(2)-1,   1-4*xi(3), ...
        4*xi(1),   4*xi(3)-4*xi(2), -4*xi(1)];

    dxdxi=Nap*X;
    Det=det(dxdxi);

    % Global derivatives [dN/dx; dN/dy].
    dNdx=[ ...
         dxdxi(2,2), -dxdxi(1,2); ...
        -dxdxi(2,1),  dxdxi(1,1)]/Det*Nap;

    eldf=12;
    inx=(2:2:eldf)'-1;
    iny=inx+1;

    B=zeros(3,eldf);
    B(1,inx)=dNdx(1,:);
    B(2,iny)=dNdx(2,:);
    B(3,inx)=dNdx(2,:);
    B(3,iny)=dNdx(1,:);
end
