function convertKnossosNmlToHoc3(a,b,c,lc,uc,filename)
fid=fopen([filename '.hoc'], 'w+');
for i=1:size(a,1)
    if c(i) > lc && c(i) < uc
        fprintf(fid,'\n{create adhoc%i}\n',i);
        fprintf(fid,'{access adhoc%i}\n',i);
        fprintf(fid,'{nseg = 1}\n');
        fprintf(fid,'{strdef color color = "White"}\n');
        fprintf(fid,'{pt3dclear()}\n');
        fprintf(fid,'{pt3dadd(%f,%f,%f,%f)}\n',b(a(i,1),1),b(a(i,1),2),b(a(i,1),3),100);
        fprintf(fid,'{pt3dadd(%f,%f,%f,%f)}\n',b(a(i,2),1),b(a(i,2),2),b(a(i,2),3),100);
    end
end
fclose(fid);
end