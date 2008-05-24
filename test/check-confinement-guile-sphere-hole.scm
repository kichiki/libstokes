; $Id: check-confinement-guile-sphere-hole.scm,v 1.1 2008/05/24 05:55:18 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "sphere+hole"
  10.0 ;; radius of the cavity at (0, 0, 0)
  1.0  ;; radius of the hole at (0, 0, 1) direction
))
