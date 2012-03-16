PRO create_aprvall, outfile=outfile, fromapvisit=fromapvisit
redux= 'v0.91'
if ~keyword_set(outfile) then outfile='$DATADIR/bovy/apogee/apall-'+redux+'.fits'
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
            ;;Add plPlugMap structure's relevant entries
            plPlugMap= mrdfits(mjdDirs[jj]+'/apPlate-b-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits',11,/silent)
            match, apvisit.fiberid, plPlugMap.fiberid, suba, subb, /sort
            newtag= {apogee_target1:0L,apogee_target2:0L,objtype:''}
            newtag= replicate(newtag,n_elements(apvisit))
            newtag[suba].apogee_target1= plPlugMap[subb].primtarget
            newtag[suba].apogee_target2= plPlugMap[subb].sectarget
            newtag[suba].objtype= plPlugMap[subb].objtype
            apvisit= struct_combine(apvisit,newtag)
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
;;flag duplicates
apogee_rem_dups, outStr, out, /flag, /structs
mwrfits, out, outfile, /create
;;add AK w/ python
cmd= "/usr/local/epd/bin/python add_ak.py "+outfile+" "+outfile
spawn, cmd
END
