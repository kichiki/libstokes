; $Id: check-ev-LJ-guile.scm,v 1.1 2008/05/24 06:07:48 kichiki Exp $
(define ev-LJ '(
 (; LJ type 1
  10.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
  (    ; 3) list of particles
   0 1 2
  )
 )
 (; LJ type 2
  8.0  ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
  2.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
  (    ; 3) list of particles
   3 4
  )
 )
))
