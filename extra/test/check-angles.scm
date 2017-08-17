; $Id: check-angles.scm,v 1.2 2008/04/17 04:20:39 kichiki Exp $
(define angles '(
  (; angle type 1
   10.0    ; 1) constant (k^{angle})
   180.0   ; 2) angle in degree (theta_0)
   0       ; 3) scale flag (0 == scaled)
           ;    in this case, the above value for k is just used.
   ((0 1 2); 3) list of triplets
    (1 2 3)
    (2 3 4)
   )
  )
  (; angle type 2
   20.0    ; 1) constant (k^{angle})
   90.0    ; 2) angle in degree (theta_0)
   1       ; 3) scale flag (1 == not scaled yet)
           ;    in this case, the potential is given by 
           ;    (k/2) * kT * (theta - theta_0)^2
   ((3 4 5); 3) list of triplets
    (4 5 6)
   )
  )
))
