function dump_ua_mat_records(outFile)
%DUMP_UA_MAT_RECORDS  Raw, untransformed JSON dump of the MAT-only record used by Fig. 9(b-d) (AAM-UA7 verifier helper).
% The independent Python verifier (verify_ua_figure_data.py) cannot read MATLAB table objects, so this helper writes
% the stored fields of paper/data/increment_plot_records.mat verbatim (no scaling, sorting, selection or derived
% quantity) together with the record SHA-256. All comparisons and transformations are recomputed in Python.
here=fileparts(mfilename('fullpath'));paper=fileparts(fileparts(here));f=fullfile(paper,'data','increment_plot_records.mat');
z=load(f,'Increment');I=z.Increment;
fid=fopen(f,'r');b=fread(fid,Inf,'*uint8');fclose(fid);md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
S=struct('record','paper/data/increment_plot_records.mat','sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[])), ...
 'has1mm',I.has1mm,'maxCommonCrackLength_mm',I.maxCommonCrackLength_mm,'metrics',table2struct(I.metrics,'ToScalar',true), ...
 'native',table2struct(I.native,'ToScalar',true),'metricsColumns',{I.metrics.Properties.VariableNames},'nativeColumns',{I.native.Properties.VariableNames});
txt=jsonencode(S,'PrettyPrint',true,'ConvertInfAndNaN',true);
% jsonencode prints 15-16 significant digits; append a full-precision copy of every numeric column (%.17g strings)
F=struct();for n=[I.metrics.Properties.VariableNames,I.native.Properties.VariableNames]
 if any(strcmp(n{1},I.metrics.Properties.VariableNames)),v=I.metrics.(n{1});else,v=I.native.(n{1});end
 if isnumeric(v),F.(n{1})=arrayfun(@(x)sprintf('%.17g',x),v(:),'UniformOutput',false);end
end
S.fullPrecision=F;txt=jsonencode(S,'PrettyPrint',true,'ConvertInfAndNaN',true);
fid=fopen(outFile,'w');fwrite(fid,txt,'char');fclose(fid);fprintf('UA7_MAT_DUMP %s\n',outFile);
end
