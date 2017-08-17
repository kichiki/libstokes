; $Id: check-confinement-guile-cylinder.scm,v 1.1 2008/05/24 05:55:52 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "cylinder"    ;; the cylinder center goes through (0, 0, 0) and (x, y, z).
  10.0          ;; radius of the cylinder
  1.0  0.0  0.0 ;; direction vector (x, y, z) of the cylinder
))
