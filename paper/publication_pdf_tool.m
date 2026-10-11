function publication_pdf_tool(args)
% Run the small PDF layout helper; PUBLICATION_PYTHON can select a runtime.
python=getenv('PUBLICATION_PYTHON');
if isempty(python)
    bundled=fullfile(getenv('USERPROFILE'),'.cache','codex-runtimes', ...
        'codex-primary-runtime','dependencies','python','python.exe');
    if isfile(bundled),python=bundled;else,python='python';end
end
script=fullfile(fileparts(mfilename('fullpath')),'normalize_retained_sensitivity_pdf.py');
parts=[{python,script},args];
for k=1:numel(parts)
    assert(~contains(parts{k},'"')&&~contains(parts{k},newline),'fig789:UnsafePath','Unsupported command path');
    parts{k}=['"' parts{k} '"'];
end
[status,message]=system(strjoin(parts,' '));
assert(status==0,'fig789:PDFToolFailure','PDF layout helper failed: %s',message);
end
