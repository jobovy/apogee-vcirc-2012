PRO create_aprvall, outfile=outfile
redux= 'v0.91'
if ~keyword_set(outfile) then outfile='$DATADIR/bovy/apogee/aprvall-'+redux+'.fits
;;Get all of the plates
basedir= '$APOGEE_ROOT/spectro/'+redux+'/plates/'
plateDirs= file_search(basedir+'*',/test_directory)
foundFirst= 0
for ii=0L, n_elements(plateDirs)-1 do begin
    mjdDirs= file_search(plateDirs[ii]+'/*',/test_directory)
    plate= (strsplit(plateDirs[ii],'/',/extract))[-1]
    for jj=0L, n_elements(mjdDirs)-1 do begin
        mjd= (strsplit(mjdDirs[jj],'/',/extract))[-1]
        filename= mjdDirs[jj]+'/apRV-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits'
        if ~file_test(filename) then continue
        if foundFirst EQ 0 then begin
            outStr= mrdfits(filename,1,/silent)
            if tag_exist(outStr,'vhelio',/quiet) then foundFirst= 1
        endif else begin
            xtra= mrdfits(filename,1,/silent)
            if tag_exist(xtra,'vhelio',/quiet) then outStr= [outStr,xtra]
        endelse
    endfor
endfor
mwrfits, outStr, outfile, /create
END
