function hires_j0826feii,redshift

   directoryname = '/Users/drupke/Box Sync/hizeahires/'
   gal = 'j0826'
   redshift = 0.603d
   readcol,directoryname+'spec/J0826_stitched_v1.txt', $
           wavelength, flux_unnorm, var, cont, flux, $
           FORMAT='D,D,D,D,D',/silent,skip=2
   error = sqrt(var)

   fittedline = 'FeII'
   reflinename = 'FeII2585'
   refloglf = 2.252d
   refmultwave = 2585.8762d
;   linenames = ['FeII2599','FeII2382','FeII2373','FeII2343']
;   loglf = [2.793d,2.882d,1.871d,2.427]
;   multwave = [2599.3959d,2382.0386d,2373.7365d,2343.4948d]
   linenames = ['FeII2599','FeII2382','FeII2373','FeII2366','FeII2343',$
                'FeII2260','FeII2249','FeII2233']
   loglf = [2.793d,2.882d,1.871d,-1.291,2.427,0.742,0.612,-1.249]
   multwave = [2599.3959d,2382.0386d,2373.7365d,2366.8674d,2343.4948d,$
               2260.0809d,2249.1795d,2233.7532d]

   bad = 1d99
   ncols = 1
   nrows = 1
   readcol, directoryname+'/fits/'+gal+'/'+gal+fittedline+'par_init.txt', $
            profileshifts, profilesig, coveringfactor, opticaldepth, $
            FORMAT='(A,D,D,D,D)',/silent
   comps=N_ELEMENTS(profileshifts)

; Finding the index to fit over
   linefitreg=[3550,4200]
   lineplotreg=[3550,4200]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   weight=1d/var
;
;   contplotreg=[3600,4200]
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
;

   relativeflux=flux
   relativeerror=error

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for multiplet fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    nabs = dblarr(ncols,nrows)+comps
    abs_zinit = dblarr(ncols,nrows,comps) + redshift
    FOR I = 0, comps-1 DO abs_zinit[*,*,I] += profileshifts[I]/refmultwave
    abs_siginit = dblarr(ncols,nrows,comps)
    FOR I = 0, N_ELEMENTS(profilesig)-1 DO $
        abs_siginit[*,*,I] += profilesig[I]
    abs_siglim = [1d,1000d]
    abs_fix = bytarr(ncols,nrows,comps,4)
;    abs_fix[0,0,*,0] = 1b
    abs_cfinit = coveringfactor
    abs_tauinit = opticaldepth

    init = {$
;
      plotindex:lineplotind,$
      fitindex:linefitind,$
      fcnfitmultiplet: 'ifsf_multipletfcn',$
      fcninitpar: 'ifsf_initmultiplet',$
      reflinename: reflinename,$
      refloglf: refloglf,$
      refmultwave: refmultwave,$
      linenames: linenames,$
      loglf: loglf,$
      multwave: multwave,$
;
      maxncomp: comps,$
      taumax: 10d,$
;
      nabs: nabs,$
      abs_cfinit: abs_cfinit,$
      abs_tauinit: abs_tauinit,$
      abs_zinit: abs_zinit,$
      abs_siginit: abs_siginit,$
      abs_siglim: abs_siglim,$
      abs_fix: abs_fix,$
;
      argspltmultiplet: {smooth: 5,$
                         xran: [[3585,3625],[3735,3765],[3770,3825],[4125,4170]]},$
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
