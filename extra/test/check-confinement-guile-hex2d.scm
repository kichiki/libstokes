; $Id: check-confinement-guile-hex2d.scm,v 1.1 2008/05/24 05:56:54 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "hex2d"
  10.0    ;; cavity radius
  1.0     ;; cylinder radius
  12.0    ;; lattice spacing
))
