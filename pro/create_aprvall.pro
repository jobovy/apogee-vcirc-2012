PRO create_aprvall, outfile=outfile, fromapvisit=fromapvisit
redux= 'v0.91'
if ~keyword_set(outfile) then outfile='$DATADIR/bovy/apogee/aprvall-'+redux+'.fits'
;;Get all of the plates
basedir= '$APOGEE_ROOT/spectro/'+redux+'/plates/'
plateDirs= file_search(basedir+'*',/test_directory)
foundFirst= 0
for ii=0L, n_elements(plateDirs)-1 do begin
    mjdDirs= file_search(plateDirs[ii]+'/*',/test_directory)
    plate= (strsplit(plateDirs[ii],'/',/extract))[-1]
    for jj=0L, n_elements(mjdDirs)-1 do begin
        mjd= (strsplit(mjdDirs[jj],'/',/extract))[-1]
        if keyword_set(fromapvisit) then begin
            apVisits= file_search(mjdDirs[jj]+'/apVisit-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'-*.fits',/test_regular)
            napvisits= n_elements(apVisits)
            if ~file_test(apVisits[0]) then continue
            apvisit= read_apvisit_header(apVisits[0])
            for kk=1L, n_elements(apVisits)-1 do begin
                apvisit= [apvisit,read_apvisit_header(apVisits[kk])]
            endfor
            if foundFirst EQ 0 then begin
                outStr= apvisit
                if tag_exist(outStr,'vhelio',/quiet) then foundFirst= 1
            endif else begin
                xtra= apvisit
                if tag_exist(xtra,'vhelio',/quiet) then outStr= [outStr,xtra]
            endelse
        endif else begin            
            filename= mjdDirs[jj]+'/apRV-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits'
            if ~file_test(filename) then continue
            if foundFirst EQ 0 then begin
                outStr= mrdfits(filename,1,/silent)
                if tag_exist(outStr,'vhelio',/quiet) then foundFirst= 1
            endif else begin
                xtra= mrdfits(filename,1,/silent)
                if tag_exist(xtra,'vhelio',/quiet) then outStr= [outStr,xtra]
            endelse
        endelse
    endfor
endfor
mwrfits, outStr, outfile, /create
END
