; $Id: check-ev-dh-guile.scm,v 1.3 2008/11/01 05:56:05 kichiki Exp $
(define ev-dh '(
  ; system parameters
  1.0e-6   ; 1) epsilon for the cut-off distance for EV_DH interaction
  298.0    ; 2) temperature [K]
  80.0     ; 3) dielectric constant of the solution
  3.07     ; 4) Debye length [nm]
  1        ; 5) flag_grid (0 == particle-particle loop, 1 == grid loop)
  (        ; 6) list of DH types
   (; DH type 1
    2.43    ; 1) nu [e/nm]
    5.00    ; 2) l0 [nm]
    (0 1 2) ; 3) list of particles
   )
   (; DH type 2
    2.00    ; 1) nu [e/nm]
    4.00    ; 2) l0 [nm]
    (3 4)   ; 3) list of particles
   )
  )
))
