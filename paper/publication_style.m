function P=publication_style()
% Physical export sizes match inline LaTeX widths: no hidden font shrinking.
P.textWidth_cm=16.2; % A4 minus the manuscript's two 24-mm margins
P.fullWidth_cm=.94*P.textWidth_cm;
P.halfWidth_cm=.485*P.textWidth_cm;
P.thirdWidth_cm=.315*P.textWidth_cm;
P.label_pt=11;P.tick_pt=9.5;P.legend_pt=11;P.annotation_pt=11;
P.lineWidth_pt=1;P.markerSize_pt=3;
P.fontName='Computer Modern'; % MATLAB LaTeX; compatible with manuscript LM
end
