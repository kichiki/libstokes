; $Id: check-confinement-guile-porous.scm,v 1.1 2008/10/31 05:50:32 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "porous"
  10.0    ;; particle radius
  20.0    ;; lattice spacing in x (2R for touching case)
))
