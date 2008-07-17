; $Id: check-excluded-volume-guile.scm,v 1.1 2008/07/17 21:47:40 kichiki Exp $
(define ev '(
 5.0     ; max distance [nm] (or in the same dimension of "length")
 ( ; for the EV 1
  0.0012 ; v [nm^3] (or in the same dimension of "length")
  0      ; fene
  1.0    ; p1 = A^{sp}, scaled spring const
  2.1    ; p2 = L_{s} / length, scaled max extension
  (0 1 2); list of particles belongs to the EV parameters
 )
 ( ; for the EV 2
  0.002  ; v [nm^3] (or in the same dimension of "length")
  1      ; fene
  19.8   ; p1 = N_{K,s}, the Kuhn steps for a spring
  106.0  ; p2 = b_{K} [nm], the Kuhn length
  (3 4)  ; list of particles belongs to the EV parameters
 )
))
