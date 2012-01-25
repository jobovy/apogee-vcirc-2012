PRO APOGEE_REM_DUPS, infile, outfile, flag=flag, structs=structs
;+
;   NAME:
;      apogee_rem_dups
;   PURPOSE:
;      remove (or flag) duplicates
;   INPUT:
;      infile - file that has the data with duplicates to be removed
;      outfile - file with the duplicates removed (or flagged)
;   KEYWORDS:
;      flag= if set, just flag the duplicates, don't remove them
;      structs= if set, input is a struct and output is a struct,
;               rather than a file
;   OUTPUT:
;      in outfile
;   HISTORY:
;      11-11 - Started - Bovy (IAS)
;-
if ~keyword_set(structs) then in= mrdfits(infile,1) else in= infile
;;first sort on snr to resolve multiple matches
if tag_exist(in,'sna') then begin
    sortindx= reverse(sort(in.sna))
endif else if tag_exist(in,'snr') then begin
    sortindx= reverse(sort(in.snr))
endif else if tag_exist(in,'vraderr') then begin
    sortVal= in.vraderr/abs(in.vrad)
    sortVal[where(sortVal EQ 0.)]= max(sortVal) ;;remove bad measurements
    sortindx= sort(sortVal)
endif
in= in[sortindx]
spherematch, in.ra, in.dec, in.ra, in.dec, 0.5/3600., m1, m2, maxmatch=0
;;first get uniqs in m1
sortindx= reverse(sort(m1))
m1= m1[sortindx]
m2= m2[sortindx]
uniqindx= uniq(m1)
m1_uniq1= m1[uniqindx]
m2_uniq1= m2[uniqindx]
;;then get uniq in m2
sortindx2= sort(m2_uniq1)
m1_uniq1= m1_uniq1[sortindx2]
m2_uniq1= m2_uniq1[sortindx2]
uniqindx2= uniq(m2_uniq1)
m1_uniq2= m1_uniq1[uniqindx2]
m2_uniq2= m2_uniq1[uniqindx2]
if keyword_set(flag) then begin
    extra=  {specprimary:0B,uniqid:0LL,specid:0LL,nvisits:0L}
    extra= replicate(extra,n_elements(in))
    extra.specid= specid(in.plate,in.mjd5,in.fiberid)
    extra[m1_uniq2].specprimary= 1B
    extra[m1_uniq1].uniqid= extra[m2_uniq1].specid
    ;;at this stage, we can record nvisits
    nvisits= lindgen(n_elements(m1_uniq1))+1L
    nvisits= nvisits[uniqindx2];;last one of series gives cumulative number of nvisits
    nvzero= nvisits[0] ;;record before shifting array
    nvisits= nvisits-shift(nvisits,1) ;;get differential number of visits
    nvisits[0]= nvzero ;;restore zero-th element
    extra[m1_uniq2].nvisits= nvisits
    extra[m1_uniq1].nvisits= extra[m2_uniq1].nvisits
    out= struct_combine(in,extra)
endif else begin
    extra=  {specprimary:1B,uniqid:0LL,specid:0LL}
    out= in[m1_uniq2]
    extra= replicate(extra,n_elements(out))
    extra.specid= specid(out.plate,out.mjd5,out.fiberid)
    extra.uniqid= extra.specid
    out= struct_combine(out,extra)
endelse
if keyword_set(structs) then outfile= out else $
  mwrfits, out, outfile, /create
END
