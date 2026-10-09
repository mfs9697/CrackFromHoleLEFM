function hash=publication_sha256(file)
% Byte identity for publication inputs; no numerical calculation.
fid=fopen(file,'rb');assert(fid>=0,'publication:MissingFile','Missing %s.',file);
cleanup=onCleanup(@()fclose(fid));
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),bytes=fread(fid,1048576,'*uint8');digest.update(bytes);end
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
