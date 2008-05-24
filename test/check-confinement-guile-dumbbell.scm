; $Id: check-confinement-guile-dumbbell.scm,v 1.1 2008/05/24 05:56:28 kichiki Exp $
(define confinement '(
  10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
  1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
  "dumbbell" ;; the origin is at the center of the cylinder
  10.0       ;; left cavity radius centered at (center1, 0, 0)
  10.0       ;; right cavity radius centered at (center2, 0, 0)
  2.0        ;; length of the cylinder
  1.0        ;; cylinder radius
))
