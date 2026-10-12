function publication_finalize_axis(ax,pdf)
% One vector-crop/PNG route for all live MATLAB panel exports.
fig=ancestor(ax,'figure');ax.Units='centimeters';
meta=struct('axes_mm',ax.Position*10,'xlim',ax.XLim,'ylim',ax.YLim, ...
    'legendPolicy','Inside own axes; collision-checked','legendBounds_mm',[]);
for lg=findall(fig,'Type','legend').'
    lg.Units='centimeters';meta.legendBounds_mm(end+1,:)=lg.Position*10;
end
file=strrep(pdf,'.pdf','.layout.json');
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(meta));fclose(fid);
publication_pdf_tool({'--crop',pdf});
end
