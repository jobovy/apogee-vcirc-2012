FUNCTION READ_APVISIT_HEADER, filename
tags= ['plate','mjd5','fiberid','location','ncombine','exptime',$
       'objid','objtype',$
       'jd','glon','glat','ra','dec','jmag','hmag','kmag','bc',$
       'vrad','vraderr','vscatter',$
       'vhelio','feh','teff','logg','xshift','ccpeak','ccpfwhm',$
       'chisq']
ntags= n_elements(tags)
format= strarr(ntags)+'F'
format[0]= 'L'
format[1]= 'L'
format[2]= 'L'
format[3]= 'L'
format[4]= 'L'
format[6]= 'A'
format[7]= 'A'
values= strarr(ntags)
out= {plate:0L}
values[*]= '0D'
indx= where(strcmp(format,'A'),cnt)
if cnt gt 0 then values[indx]= "''"
indx= where(strcmp(format,'I'),cnt)
if cnt gt 0 then values[indx]= '0'
indx= where(strcmp(format,'L'),cnt)
if cnt gt 0 then values[indx]= '0L'
out= struct_addtags(out,tags[1:-1],values[1:-1])
tnames=TAG_NAMES(out)
;;read
hdr= headfits(filename)
for ii=0L, ntags-1 do begin
    tindex=WHERE(STRCMP(tnames,tags[ii],/fold_case) EQ 1)
    out.(tindex)= sxpar(hdr,tags[ii])
endfor
return, out
END
