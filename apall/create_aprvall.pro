PRO create_aprvall, outfile=outfile, redux=redux, $
                    silent=silent, aspcap=aspcap, abredux=abredux
if ~keyword_set(redux) then redux= 'v1'
if ~keyword_set(abredux) then abredux= 'v0.3_1'
if ~keyword_set(outfile) then begin
    if keyword_set(aspcap) then $
      outfile='$DATADIR/bovy/apogee/apall-1d-'+redux+'-aspcap-'+abredux+'.fits' $
    else $
      outfile='$DATADIR/bovy/apogee/apall-'+redux+'.fits'
endif
;;Get all of the plates
;;OLD basedir= '$APOGEE_ROOT/spectro/'+redux+'/plates/'
basedir= '$APOGEE_ROOT/spectro/redux/'+redux+'/plates/'
aspcapdir= '/mount/hydra4/jb2777/apogee_data/aspcap/'+abredux+'/'
plateDirs= file_search(basedir+'*',/test_directory)
foundFirst= 0
for ii=0L, n_elements(plateDirs)-1 do begin
    mjdDirs= file_search(plateDirs[ii]+'/*',/test_directory)
    plate= (strsplit(plateDirs[ii],'/',/extract))[-1]
    if ~keyword_set(silent) then print, format = '("Working on plate ",i4,a1,$)' $
      , plate, string(13b)
    for jj=0L, n_elements(mjdDirs)-1 do begin
        mjd= (strsplit(mjdDirs[jj],'/',/extract))[-1]
        apVisits= file_search(mjdDirs[jj]+'/apVisit-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'-*.fits',/test_regular)
        napvisits= n_elements(apVisits)
        if ~file_test(apVisits[0]) then continue
        apvisit= read_apvisit_header(apVisits[0])
        for kk=1L, n_elements(apVisits)-1 do begin
            apvisit= [apvisit,read_apvisit_header(apVisits[kk])]
        endfor
        ;;Add plPlugMap structure's relevant entries
        plPlugMap= mrdfits(mjdDirs[jj]+'/apPlate-b-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits',11,/silent,status=status)
        if status lt 0 then begin
            print, "plPugMap for this version not found, trying v0.91 ..."
            plPlugMap= mrdfits(repstr(mjdDirs[jj]+'/apPlate-b-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits','v1','v0.91'),11,/silent,status=status)
            if status lt 0 then print, "plPlugMap not found ..."
        endif
        match, apvisit.fiberid, plPlugMap.fiberid, suba, subb, /sort
        newtag= {apogee_target1:0L,apogee_target2:0L,objtype:''}
        newtag= replicate(newtag,n_elements(apvisit))
        newtag[suba].apogee_target1= plPlugMap[subb].primtarget
        newtag[suba].apogee_target2= plPlugMap[subb].sectarget
        newtag[suba].objtype= plPlugMap[subb].objtype
        apvisit= struct_combine(apvisit,newtag)
        if keyword_set(aspcap) then begin
           
            ;;Add aspcapPlate's relevant entries
            aspcapPlateFile= aspcapdir+'aspcapPlate-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits'
            newtag= {param:dblarr(7),param_cova:dblarr(7,7),$
                     sn:-9999.00D,sn_median:-9999.00D,$
                     aspcap_feh:-9999.00D,afe:-9999.00D,$
                     cfe:-9999.00D,nfe:-9999.00D,$
                     aspcap_logg:-9999.00D,aspcap_teff:-9999.00D}
            newtag= replicate(newtag,n_elements(apvisit))
            if ~file_test(aspcapPlateFile) then begin
                print, "ASPCAP for plate "+strtrim(string(plate),2)+", MJD "+strtrim(string(mjd),2)+" not available"
            endif else begin
                aspcapPlate= mrdfits(aspcapPlateFile,1,/silent)
                ;;match
                match, apvisit.fiberid, aspcapPlate.fiberid, suba, subb, /sort
                newtag[suba].param= aspcapPlate[subb].param
                newtag[suba].param_cova= aspcapPlate[subb].param_cova
                newtag[suba].sn= aspcapPlate[subb].sn
                newtag[suba].sn_median= aspcapPlate[subb].sn_median
                newtag[suba].aspcap_feh= aspcapPlate[subb].param[3]
                newtag[suba].aspcap_teff= aspcapPlate[subb].param[0]
                newtag[suba].aspcap_logg= aspcapPlate[subb].param[1]
                newtag[suba].cfe= aspcapPlate[subb].param[4]
                newtag[suba].nfe= aspcapPlate[subb].param[5]
                newtag[suba].afe= aspcapPlate[subb].param[6]
            endelse
            apvisit= struct_combine(apvisit,newtag)
        endif
        if foundFirst EQ 0 then begin
            outStr= apvisit
            if tag_exist(outStr,'vhelio',/quiet) then foundFirst= 1
        endif else begin
            xtra= apvisit
            if tag_exist(xtra,'vhelio',/quiet) then outStr= [outStr,xtra]
        endelse
    endfor
endfor
;;rename reduction Feh etc.
if keyword_set(aspcap) then begin
    outStr= rename_tags(outStr,['feh','logg','teff'],$
                        ['feh1D','logg1D','teff1D'])
    outStr= rename_tags(outStr,['aspcap_feh','aspcap_logg','aspcap_teff'],$
                        ['feh','logg','teff'])
endif
;;flag duplicates
apogee_rem_dups, outStr, out, /flag, /structs
mwrfits, out, outfile, /create
;;add AK w/ python
cmd= "/usr/local/epd/bin/python add_ak.py "+expand_path(outfile)+" "+expand_path(outfile)
spawn, cmd
;;add averages w/ python
cmd= "/usr/local/epd/bin/python add_avg_ab.py "+expand_path(outfile)+" "+expand_path(outfile)
spawn, cmd
END
