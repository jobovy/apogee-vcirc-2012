PRO APOGEE_REM_DUPS, infile, outfile
in= mrdfits(infile,1)
;;first sort on snr to resolve multiple matches
if tag_exist(in,'sna') then begin
    sortindx= reverse(sort(in.sna))
endif else if tag_exist(in,'snr') then begin
    sortindx= reverse(sort(in.snr))
endif else if tag_exist(in,'vraderr') then begin
    sortindx= sort(in.vraderr)
endif
in= in[sortindx]
spherematch, in.ra, in.dec, in.ra, in.dec, 0.5/3600., m1, m2, maxmatch=0
;;first get uniqs in m1
sortindx= sort(m1)
m1= m1[sortindx]
m2= m2[sortindx]
uniqindx= uniq(m1)
m1= m1[uniqindx]
m2= m2[uniqindx]
;;then get uniq in m2
sortindx2= sort(m2)
m1= m1[sortindx2]
m2= m2[sortindx2]
uniqindx2= uniq(m2)
m1= m1[uniqindx2]
m2= m2[uniqindx2]
out= in[m1]
mwrfits, out, outfile, /create
END
