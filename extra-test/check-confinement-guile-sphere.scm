; $Id: check-confinement-guile-sphere.scm,v 1.1 2008/05/24 05:54:10 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "sphere"
  10.0 ;; radius of the cavity at (0, 0, 0)
))
