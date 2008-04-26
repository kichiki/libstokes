; $Id: check-ev-dh-guile.scm,v 1.1 2008/04/26 05:05:54 kichiki Exp $
(define ev-dh '(
  ; system parameters
  4.0      ; 1) max distance for EV_DH interaction [nm]
  298.0    ; 2) temperature [K]
  80.0     ; 3) dielectric constant of the solution
  3.07     ; 4) Debye length [nm]
  (        ; 5) list of chain types
   (; chain type 1
    2.43    ; 1) nu [e/nm]
    5.00    ; 2) l0 [nm]
    (0 1 2) ; 3) list of particles
   )
   (; chain type 2
    2.00    ; 1) nu [e/nm]
    4.00    ; 2) l0 [nm]
    (3 4)   ; 3) list of particles
   )
  )
))
