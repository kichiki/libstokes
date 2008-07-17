; $Id: check-bead-rod-guile.scm,v 1.1 2008/07/17 03:05:32 kichiki Exp $
(define constraints '(
 ; system parameters
 1.0e-6    ; 1) tolerance
 "nitsol"  ; 2) scheme for solving nonlinear equations
           ;    "linear" for iterative scheme in linear approximation
           ;    "nitsol" for Newton-GMRES scheme by NITSOL library
 ; the following is for each constraint
 (         ; 3) constraint type 1
  5.0      ; 3-1) distance [nm]
  (        ; 3-2) list of particle-pairs
   (0 1)
   (1 2)
   (2 3)
 ))
 (         ; 4) constraint type 2
  10.0     ; 4-1) distance [nm]
  (        ; 4-2) list of particle-pairs
   (3 4)
   (4 5)
 ))
))
