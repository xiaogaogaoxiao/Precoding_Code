 
function Save_link=Save_link_file(readFile,writeFile)
fid = fopen(readFile);
fid_write = fopen(writeFile,'w');
while 1
     tline = fgetl(fid);
    if ~ischar(tline),
        break, 
    end
    %tline = fgetl(fid);
    fprintf(fid_write,'%s\n',tline);
    %disp(tline)
end
fclose(fid);
fclose(fid_write);
Save_link=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 