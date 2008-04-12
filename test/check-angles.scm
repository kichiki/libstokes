; $Id: check-angles.scm,v 1.1 2008/04/12 19:17:23 kichiki Exp $
(define angles '(
  (; angle type 1
   10.0    ; 1) constant (k^{angle})
   180.0   ; 2) angle in degree (theta_0)
   ((0 1 2); 3) list of triplets
    (1 2 3)
    (2 3 4)
   )
  )
  (; angle type 2
   20.0    ; 1) constant (k^{angle})
   90.0    ; 2) angle in degree (theta_0)
   ((3 4 5); 3) list of triplets
    (4 5 6)
   )
  )
))
