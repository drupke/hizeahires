function hires_j0826mgii,redshift,rootdir

   directoryname = '/Users/drupke/Box Sync/hizeahires/'
   gal = 'j0826'
   redshift = 0.603d
   readcol,directoryname+'spec/J0826_stitched_v1.txt', $
           wavelength_air, flux_unnorm, var, cont, flux, $
           FORMAT='D,D,D,D,D',/silent,skip=2
;  normalize variance
   error = sqrt(var)/(flux_unnorm/flux)
;  convert to vacuum wavelengths
   airtovac, wavelength_air, wavelength

   fittedline = 'MgII'
   restred = 2803.530d ; vacuum
   
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

; Finding the index to fit over
   linefitreg=[4450,4550]
   lineplotreg=[4450,4510]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   weight=1d/error^2d

;   contplotreg=[4400,4600]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;      VALUE_LOCATE(wavelength,contplotreg[1])]
;   set_plot,'z'
;   cgplot, wavelength, flux, XRAN=contplotreg, $
;           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
;           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
;           xtit='Wavelength ($\Angstrom$)',$
;           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;   img = cgsnapshot(filename=directoryname+'/fits/'+gal+'/'+gal+fittedline+$
;                    '_continuum',/jpeg,/nodialog,quality=100)

;   nct = contplotind[1]-contplotind[0]+1
;   nflux = dblarr(nct)
;   nerr = dblarr(nct)
;   openw,outlun,directoryname+'/'+gal+fittedline+'_contfit.txt',/get_lun
;   for i=0,nct-1 do $
;      printf,outlun,wavelength[contplotind[0]+i],flux[contplotind[0]+i],$
;             error[contplotind[0]+i],$
;             flux[contplotind[0]+i]/continuum[contplotind[0]+i],$
;             error[contplotind[0]+i]/continuum[contplotind[0]+i],$
;             format='(D12.4,E12.4,E12.4,D8.4,D8.4)'
;   free_lun,outlun

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
;    doubletabs_fix[0,0,*,0] = 1b
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
    doubletem_fix[0,0,*,0] = 1b
    ;doubletem_fix[0,0,*,2] = 1b
    doubletem_fix[0,0,*,3] = 1b
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
;
      ndoubletabs: ndoubletabs,$
      doubletabs_cfinit: doubletabs_cfinit,$
      doubletabs_tauinit: doubletabs_tauinit,$
      doubletabs_zinit: doubletabs_zinit,$
      doubletabs_siginit: doubletabs_siginit,$
      doubletabs_siglim: doubletabs_siglim,$
      doubletabs_fix: doubletabs_fix,$
;
      taumax: 10d,$
      mcniter: 1000,$
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
