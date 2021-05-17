function hires_j1341mgii,redshift

   directoryname = '/Users/drupke/Box Sync/hizeahires/'
   gal = 'j1341'
   redshift = 0.661d
   readcol,directoryname+'spec/J1341_stitched_v1.txt', $
           wavelength, flux_unnorm, var, cont, flux, $
           FORMAT='D,D,D,D,D',/silent,skip=2
   error = sqrt(var)

   fittedline = 'MgII'
   restred = 2802.704d
   
   bad = 1d99
   ncols = 1
   nrows = 1
   readcol, directoryname+'/fits/'+gal+'/'+gal+fittedline+'par_init.txt', $
            profileshifts, profilesig, coveringfactor, opticaldepth, $
            FORMAT='(A,D,D,D,D)',/silent
   comps=N_ELEMENTS(profileshifts)
   readcol, directoryname+'/fits/'+gal+'/'+gal+fittedline+'par_init_em.txt', $
            profileshifts_em, profilesig_em, emflux, emratio, $
            FORMAT='(A,D,D,D,D)',/silent
   comps_em=N_ELEMENTS(profileshifts_em)
;   comps_em = 0

; Finding the index to fit over
   linefitreg=[4605,4675]
   lineplotreg=[4605,4675]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   
   weight=1d/var
   relativeflux=flux
   relativeerror=error

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for doublet fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + redshift
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/restred
    doubletabs_siginit = dblarr(ncols,nrows,comps)
    FOR I = 0, N_ELEMENTS(profilesig)-1 DO $
        doubletabs_siginit[*,*,I] += profilesig[I]
    doubletabs_siglim = [1d,1000d]
    doubletabs_fix = bytarr(ncols,nrows,comps,4)
;    doubletabs_fix[0,0,*,0]=1b
    doubletabs_cfinit = coveringfactor
    doubletabs_tauinit = opticaldepth

    ndoubletem = dblarr(ncols,nrows)+comps_em
    doubletem_zinit = dblarr(ncols,nrows,comps_em) + redshift
    FOR I = 0, comps_em-1 DO doubletem_zinit[*,*,I] += profileshifts_em[I]/restred
    doubletem_siginit = dblarr(ncols,nrows,comps_em)
    FOR I = 0, N_ELEMENTS(profilesig_em)-1 DO $
        doubletem_siginit[*,*,I] += profilesig_em[I]
    doubletem_siglim = [1d,1000d]
    doubletem_fix = bytarr(ncols,nrows,comps_em,4)
    doubletem_fix[*,*,*,1:2] = 1b
    doubletem_finit = emflux
    doubletem_rinit = emratio

    init = {$
;
      plotindex:lineplotind,$
      fitindex:linefitind,$
      fcnfitdoublet: 'ifsf_doubletfcn',$
      fcninitpar: 'ifsf_initdoublet',$
;
      maxncomp: comps,$
      taumax: 10,$
;
      ndoubletabs: ndoubletabs,$
      doubletabs_cfinit: doubletabs_cfinit,$
      doubletabs_tauinit: doubletabs_tauinit,$
      doubletabs_zinit: doubletabs_zinit,$
      doubletabs_siginit: doubletabs_siginit,$
      doubletabs_siglim: doubletabs_siglim,$
      doubletabs_fix: doubletabs_fix,$
;
      ndoubletem: ndoubletem,$
      doubletem_zinit: doubletem_zinit,$
      doubletem_siginit: doubletem_siginit,$
      doubletem_finit: doubletem_finit,$
      doubletem_rinit: doubletem_rinit,$
      doubletem_siglim: doubletem_siglim,$
      doubletem_fix: doubletem_fix,$
;
      galaxy: gal, $
      zsys_gas: redshift,$
      outdir: directoryname+'/fits/'+gal+'/',$
      wavelength: wavelength, $
      relativeflux: relativeflux, $
      error: relativeerror, $
      continuum: flux/flux, $
      flux: flux $
    }

  return,init

END
