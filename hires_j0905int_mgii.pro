function hires_j0905int_mgii,redshift

   directoryname = '/Users/drupke/Box Sync/hizeahires/'
   gal = 'j0905int'
   redshift = 0.412d
   readcol,directoryname+'spec/J0905_stitched_v1.txt', $
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

; Finding the index to fit over
   linefitreg=[3940,3970]
   lineplotreg=[3940,3970]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   weight=1d/var

;   contplotreg=lineplotreg
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;      VALUE_LOCATE(wavelength,contplotreg[1])]
;  set_plot,'z'
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


    ndoubletem = 0

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
;
      ndoubletem: ndoubletem,$
;
      taumax: 10d,$
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
