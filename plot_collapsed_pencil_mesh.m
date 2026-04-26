function plot_collapsed_pencil_mesh(Mc, varargin)
%PLOT_COLLAPSED_PENCIL_MESH
% Plot the collapsed mesh and highlight the collapsed crack faces.
%
% Usage:
%   plot_collapsed_pencil_mesh(Mc)
%   plot_collapsed_pencil_mesh(Mc, 'ShowOriginalFaces', true)
%
% Input:
%   Mc  struct from collapse_pencil_faces_to_midline
%
% Name-value options:
%   'ShowOriginalFaces' : logical (default false)
%   'ShowNodeNumbers'   : logical (default false)

    ip = inputParser;
    addParameter(ip, 'ShowOriginalFaces', false, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'ShowNodeNumbers', false, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});

    showOriginal = ip.Results.ShowOriginalFaces;
    showNums     = ip.Results.ShowNodeNumbers;

    must(Mc, 'p');
    must(Mc, 't');
    must(Mc, 'crack');
    must(Mc.crack, 'Pmid');
    must(Mc.crack, 'upperNodes');
    must(Mc.crack, 'lowerNodes');
    must(Mc.crack, 'tipNode');

    p  = Mc.p;
    t  = Mc.t;
    Pm = Mc.crack.Pmid;

    nU = Mc.crack.upperNodes;
    nL = Mc.crack.lowerNodes;
    nT = Mc.crack.tipNode;

    figure('Name', 'Collapsed pencil mesh', 'Color', 'w'); clf
    hold on; axis equal; box on

    % Mesh
    triplot(t(:,1:3), p(:,1), p(:,2), 'Color', [0.78 0.78 0.78]);

    % Midline
    plot(Pm(:,1), Pm(:,2), 'k-', 'LineWidth', 1.8);

    % Collapsed crack faces
    plot(p(nU,1), p(nU,2), 'ro-', 'LineWidth', 1.2, 'MarkerSize', 4);
    plot(p(nL,1), p(nL,2), 'bo-', 'LineWidth', 1.2, 'MarkerSize', 4);

    % Tip node
    plot(p(nT,1), p(nT,2), 'kp', 'MarkerSize', 10, 'LineWidth', 1.5);

    % Optional overlay of original face locations
    if showOriginal && isfield(Mc, 'p0') && ~isempty(Mc.p0)
        p0 = Mc.p0;
        plot(p0(nU,1), p0(nU,2), 'r--', 'LineWidth', 0.8);
        plot(p0(nL,1), p0(nL,2), 'b--', 'LineWidth', 0.8);
    end

    % Optional node numbers
    if showNums
        for k = 1:numel(nU)
            text(p(nU(k),1), p(nU(k),2), sprintf('U%d', nU(k)), ...
                'Color', [0.7 0 0], 'FontSize', 8, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        end
        for k = 1:numel(nL)
            text(p(nL(k),1), p(nL(k),2), sprintf('L%d', nL(k)), ...
                'Color', [0 0 0.7], 'FontSize', 8, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
        text(p(nT,1), p(nT,2), sprintf('T%d', nT), ...
            'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end

    xlabel('x')
    ylabel('y')
    title('Collapsed T3 mesh with crack faces mapped to the midline')
end


function must(S, field)
    if ~isfield(S, field) || isempty(S.(field))
        error('plot_collapsed_pencil_mesh:MissingField', ...
            'Required field "%s" is missing or empty.', field);
    end
end