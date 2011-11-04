;;Plot aprv values versus apvisit values
PRO aprv_vs_apvisit, plate, mjd, tag=tag, plotfilename=plotfilename
if ~keyword_set(tag) then tag= 'vhelio'
redux= 'v0.91'
;;Get plate/MJD
basedir= '$APOGEE_ROOT/spectro/'+redux+'/plates/'
platedir= basedir+strtrim(string(plate,format='(I4)'),2)+'/'
mjddir= platedir+strtrim(string(mjd,format='(I5)'),2)
;;check whether there is an apRV file
filename= mjddir+'/apRV-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'.fits'
if ~file_test(filename) then begin
    print, "No apRV file found ..."
    print, "Returning ..."
    return
endif
;;read apRV
aprv= mrdfits(filename,1)
;;read all apVisit files
apVisits= file_search(mjddir+'/apVisit-'+strtrim(string(plate,format='(I4)'),2)+'-'+strtrim(string(mjd,format='(I5)'),2)+'-*.fits',/test_regular)
apvisit= read_apvisit_header(apVisits[0])
for ii=1L, n_elements(apVisits)-1 do begin
    apvisit= [apvisit,read_apvisit_header(apVisits[ii])]
endfor
;;plot
tnames=TAG_NAMES(aprv)
aprv_index=WHERE(STRCMP(tnames,tag,/fold_case) EQ 1)
tnames=TAG_NAMES(apvisit)
apvisit_index=WHERE(STRCMP(tnames,tag,/fold_case) EQ 1)
if keyword_set(plotfilename) then k_print, filename=plotfilename
djs_plot, apvisit.(apvisit_index), aprv.(aprv_index)-apvisit.(apvisit_index), $
  psym=3, $
  xtitle='apVisit '+tag, $
  ytitle='apRV '+tag+' - apVisit '+tag, $
  title='plate= '+strtrim(string(plate,format='(I4)'),2)+', MJD='+strtrim(string(mjd,format='(I5)'),2)
if keyword_set(plotfilename) then k_end_print
END
