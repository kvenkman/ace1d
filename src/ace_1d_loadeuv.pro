; This code section defines quantities needed to run and use the EUVAC solar spectrum model
; This includes: 
; 1. EUVAC parameters (wavelength bins, reference spectrum, a-factors)
; 2. Absorption cross sections for various species
; 3. Branching ratios and scaling factors for various processes (dissociation, ionization, photoelectron dissociation/ionization etc.)
; Values primarily obtained from the NCAR TIE-GCM

; Wavelength bins in angstroms
solspec.wave1 = [1700.00, 1650.00, 1600.00, 1550.00, 1500.00,$
                 1450.00, 1400.00, 1350.00, 1300.00, 1250.00,$
                 1200.00, 1215.67, 1150.00, 1100.00, 1050.00,$
                 1027.00,  987.00,  975.00,  913.00,  913.00,$
                  913.00,  798.00,  798.00,  798.00,  650.00,$
                  650.00,  540.00,  320.00,  290.00,  224.00,$
                  155.00,   70.00,   32.00,   18.00,    8.00,$
                    4.00,    0.50]
                  
solspec.wave2 = [1750.00, 1700.00, 1650.00, 1600.00, 1550.00,$
                 1500.00, 1450.00, 1400.00, 1350.00, 1300.00,$
                 1250.00, 1215.67, 1200.00, 1150.00, 1100.00,$
                 1050.00, 1027.00,  987.00,  975.00,  975.00,$
                  975.00,  913.00,  913.00,  913.00,  798.00,$
                  798.00,  650.00,  540.00,  320.00,  290.00,$
                  224.00,  155.00,   70.00,   32.00,   18.00,$
                    8.00,    4.00]
;
; Reference solar spectrum, photon /cm^2 /s
;
solspec.reference = [3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,$
                     5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,$
                     2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,$
                     8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,$
                     4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,$
                     3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,$
                     1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,$
                     8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,$
                     5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,$
                     5.010e+01]
;
; scaling factor A as defined in EUVAC model
;
solspec.afac = [5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,$
                1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,$
                2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,$
                2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,$
                5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,$
                4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,$
                3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,$
                7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,$
                1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,$
                6.240e-01]

; Cross section data                
;
; Total absorption cross sections (10^-18 cm^2):
; 
; O absorption cross sections:
;

xsec.abs(0,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00,$
               8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00,$
               3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01,$
               1.70e-02, 2.27e-03]
;
; O2 absorption cross sections:
;
xsec.abs(1,*) = [$
               5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01,$
               1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01,$
               1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00,$
               1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01,$
               1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01,$
               2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01,$
               7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01,$
               3.40e-02, 4.54e-03]
;
; N2 absorption cross sections:
;
xsec.abs(2,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01,$
               2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01,$
               2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00,$
               5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01,$
               2.01e-02, 2.53e-03]
;
; CO
;
xsec.abs(3,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 1.92e+00, 3.53e+00, 5.48e+00, 8.02e+00,$
               1.00e+01, 1.17e+01, 1.10e+01, 1.25e+01, 1.25e+01,$
               1.54e+01, 1.47e+01, 1.95e+01, 2.01e+01, 2.01e+01,$
               2.01e+01, 2.18e+01, 2.18e+01, 2.18e+01, 2.23e+01,$
               2.22e+01, 2.38e+01, 1.86e+01, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
;  
; aCO2
;
xsec.abs(4,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 4.42e+00, 7.51e+00, 1.10e+01, 1.50e+01,$
               1.79e+01, 2.12e+01, 2.00e+01, 2.34e+01, 2.34e+01,$
               2.57e+01, 2.49e+01, 2.83e+01, 2.93e+01, 2.93e+01,$
               2.93e+01, 3.21e+01, 3.21e+01, 3.21e+01, 3.20e+01,$
               2.57e+01, 3.36e+01, 2.02e+01, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
;  
;O3
;
xsec.abs(5,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 1.25e+01, 9.20e+00,$
               9.20e+00, 9.20e+00, 9.16e+00, 9.50e+00, 9.50e+00,$
               9.50e+00, 1.47e+01, 1.47e+01, 1.47e+01, 2.74e+01,$
               2.02e+01, 3.33e+01, 7.74e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
;  
; bCO2
;
xsec.abs(6,*) = [$
               5.00e-02, 1.00e-01, 1.50e-01, 3.00e-01, 4.00e-01,$
               5.50e-01, 5.00e-01, 5.00e-01, 8.00e-01, 5.00e-01,$
               2.00e-01, 0.00e+00, 0.00e+00, 1.85e+01, 1.48e+01,$
               1.42e+01, 1.66e+01, 3.98e+01, 7.41e+01, 7.41e+01,$
               7.41e+01, 1.74e+01, 1.74e+01, 1.74e+01, 3.23e+01,$
               3.31e+01, 3.33e+01, 2.50e+01, 2.34e+01, 2.33e+01,$
               1.11e+01, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
;  
; H2O
;
xsec.abs(7,*) = [$
               5.00e+00, 5.00e+00, 5.00e+00, 3.00e+00, 1.50e+00,$
               8.00e-01, 8.00e-01, 1.10e+00, 5.00e+00, 8.00e+00,$
               8.00e+00, 0.00e+00, 4.44e+00, 4.44e+00, 4.44e+00,$
               1.41e+01, 2.46e+01, 1.10e+01, 1.85e+01, 1.85e+01,$
               1.85e+01, 2.36e+01, 2.36e+01, 2.36e+01, 3.98e+01,$
               2.33e+01, 2.28e+01, 2.45e+01, 2.22e+01, 1.91e+01,$
               1.11e+01, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
;  
; NO
; 
xsec.abs(8,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 1.85e+00, 2.04e+00,$
               2.04e+00, 0.00e+00, 2.41e+00, 3.70e+00, 6.48e+00,$
               3.70e+00, 7.82e+00, 1.98e+01, 2.41e+01, 2.41e+01,$
               2.41e+01, 1.55e+01, 1.55e+01, 1.55e+01, 1.36e+01,$
               1.40e+01, 2.00e+01, 2.02e+01, 2.30e+01, 2.40e+01,$
               2.22e+01, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
               
xsec.abs = xsec.abs * 1E-18 ; cm^2

; The three major species' ionization branching ratio (fraction of total absorption):
;
; O
;
xsec.br_pi(0,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00]
;
; O2 
;
xsec.br_pi(1,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01,$
               7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01,$
               6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00]
;
; N2
;
xsec.br_pi(2,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01,$
               6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00]
;
; photon ionization branching ratio for O+(2p),O+(2d),O+(4s) 
; (off O photon ionization)

; note that the sum of these 3 branching ratio is the same as xsec.br_pi[0,*]

; O+(2p)
;
branching.pio_op_2p = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               8.56e-03, 2.52e-01, 2.60e-01, 2.46e-01, 2.41e-01,$
               2.33e-01, 2.27e-01, 2.26e-01, 2.24e-01, 2.24e-01,$
               2.24e-01, 2.24e-01]
;
; O+(2d)
;
branching.pio_op_2d = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 6.98e-02,$
               3.37e-01, 4.51e-01, 4.24e-01, 4.03e-01, 4.02e-01,$
               3.92e-01, 3.77e-01, 3.74e-01, 3.78e-01, 3.78e-01,$
               3.78e-01, 3.78e-01]
;
; O+(4s)
;
branching.pio_op_4s = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 9.30e-01,$
               6.55e-01, 2.98e-01, 3.17e-01, 3.46e-01, 3.50e-01,$
               3.67e-01, 3.89e-01, 3.93e-01, 3.90e-01, 3.90e-01,$
               3.90e-01, 3.90e-01]

;
; O2 photon dissociation branching ratio
;
branching.pdo2 = [$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01,$
               2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01,$
               3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]

; this is crudely taken from Lee et al 1973 - 
; it should be done better, and needs to account for PE production of O1D?
branching.pdo2_o1dyield = [$     
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 0.50e+00, 0.50e+00, 0.50e+00,$
               0.50e+00, 0.50e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00]
               
; n2 photon dissociation branching ratio
branching.pdn2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 5.71e-01,$
               3.20e-01, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00] 
;
; O2 photon dissociative ionization branching ratio
;
branching.pdio2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               5.35e-04, 1.08e-01, 2.40e-01, 3.51e-01, 3.76e-01,$
               4.47e-01, 6.53e-01, 8.92e-01, 1.00e+00, 1.00e+00,$
               1.00e+00, 1.00e+00]
;
; n2 photon dissociative ionization branching ratio
;
branching.pdin2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 4.72e-03, 9.27e-02, 2.46e-01,$
               2.53e-01, 2.49e-01, 2.82e-01, 9.60e-01, 9.60e-01,$
               9.60e-01, 9.60e-01]
;
; n(4s) photoionization cross section
;
        xsec.sigin4s = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 3.24e+00, 2.48e+00, 2.11e+00, 1.12e+01,$
               1.18e+01, 1.16e+01, 9.35e+00, 6.43e+00, 4.80e+00,$
               2.55e+00, 6.81e-01, 1.66e-01, 5.68e-01, 7.05e-02,$
               1.00e-02, 1.27e-03]


      xsec.sigin4s = xsec.sigin4s * 1E-18

; The three major species' photoelectron ionization scaling factor
; off its photon ionization rate
;
; O
;
xsec.pepiscale(0,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 1.27e-01, 4.18e-01, 6.94e-01,$
               1.09e+00, 2.19e+00, 4.99e+00, 7.14e+01, 2.36e+01,$
               5.06e+01, 2.17e+02]
;
; O2
;
xsec.pepiscale(1,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 2.38e-02, 1.05e-01, 2.42e-01,$
               5.79e-01, 1.61e+00, 4.27e+00, 6.00e+01, 2.03e+01,$
               5.02e+01, 2.11e+02]
;
; N2
;
xsec.pepiscale(2,*) = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 3.07e-02, 1.78e-01, 3.61e-01,$
               9.33e-01, 2.86e+00, 7.79e+00, 1.08e+01, 3.22e+01,$
               8.09e+01, 3.43e+02]
;
; photoelectron impact ionization scaling factor for O+(2P),O+(2D),O+(4S)
; (off O photon ionization). These sum to pepiscale(0,*)
;
; O+(2P)
;
branching.peio_op2p = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 9.15e-03, 6.12e-02, 1.16e-01,$
               2.03e-01, 4.36e-01, 1.01e+00, 1.46e+01, 4.77e+00,$
               1.10e+01, 4.74e+01]
;
; O+(2D)      
;
branching.peio_op2d = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 3.39e-02, 1.48e-01, 2.53e-01,$
               4.18e-01, 8.53e-01, 1.96e+00, 2.82e+01, 9.36e+00,$
               2.07e+01, 8.85e+01]
;
; O+(4S)      
;
branching.peio_op4s = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 8.36e-02, 2.09e-01, 3.25e-01,$
               4.70e-01, 9.02e-01, 2.02e+00, 2.86e+01, 9.42e+00,$
               1.89e+01, 8.12e+01]
;
; photoelectron dissociative ionization scaling factor  of O2
;
xsec.pepiscale_pedio2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 9.60e-04, 1.43e-02, 5.23e-02,$
               1.63e-01, 5.21e-01, 1.44e+00, 2.03e+01, 6.98e+00,$
               1.79e+01, 7.61e+01]
;
; photoelectron dissociative ionization scaling factor of N2
;
xsec.pepiscale_pedin2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 1.84e-04, 8.49e-03, 3.66e-02,$
               1.46e-01, 5.71e-01, 1.65e+00, 2.29e+00, 6.95e+00,$
               1.83e+01, 7.87e+01]
;
; photoelectron dissociation scaling factor of N2
;
xsec.pepiscale_pedn2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 1.57e-01, 5.15e-01, 7.64e-01,$
               1.37e+00, 2.91e+00, 6.53e+00, 9.05e+00, 2.53e+01,$
               5.21e+01, 2.45e+02]
;
; photoelectron dissociation scaling factor of O2
;
xsec.pepiscale_pedo2 = [$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,$
               0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01,$
               1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01,$
               2.03e+01, 8.79e+01]
 

