function Report=redraw_all_figures(varargin)
% Reproducible plot-only vector pipeline; approved sources never overwritten.
% Issue #78 requires author inspection/approval before publication promotion.
Report=review_issue78_figures(varargin{:});
fprintf('31 candidate vector panels: %s\n',Report.figureRoot);
fprintf('Publication promotion remains pending author approval (Issue #78).\n');
end
