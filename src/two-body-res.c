/* subroutines of exact 2-body resistance scalar functions
 * Copyright (C) 1999 Kengo ICHIKI <kengo@caltech.edu>
 * $Id: two-body-res.c,v 1.1 1999/08/25 01:11:56 ichiki Exp $
 */
#include <math.h> /* log() */

#include "two-body-res.h"


/* calc scalar functions of two-body resistance
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   res [22] : scalar functions
 */
void
scalar_two_body_res (double s, double *res)
{
  static double
    one3   =  1.0/3.0, 
    two3   =  2.0/3.0, 
    one12  =  1.0/12.0, 
    two15  =  2.0/15.0, 
    one24  =  1.0/24.0, 
    one30  =  1.0/30.0, 
    g3yh   =  137.0/1500.0, 
    g5yh   =  113.0/1500.0;
  double
    x11a, x12a, y11a, y12a, 
    y11b, y12b, 
    x11c, x12c, y11c, y12c, 
    x11g, x12g, y11g, y12g, 
    y11h, y12h, 
    x11m, x12m, y11m, y12m, z11m, z12m;
  double
    s1, s2, s3, 
    s4, s5, ls1, ls2;
  
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  x11a = 
    0.25 / s1
    - 0.225 * ls1
    - (3.0 / 112.0) * s1 * ls1
    + 1.0/*XAf [0]*/
    - 0.25
    + s4 *
    (2.428571428571428e-01 + s4 *
     (2.267857142857142e-01 + s4 *
      (-1.811160714285714e+00 + s4 *
       (-4.027901785714290e-01 + s4 *
	(1.196037388392857e+01 + s4 *
	 (2.678687918526786e+01 + s4 *
	  (-1.362215775470351e+01 + s4 *
	   (-2.352605777039820e+02 + s4 *
	    (-7.253811945960639e+02 + s4 *
	     (-9.635069355537769e+02 + s4 *
	      (2.071450316300950e+03 + s4 *
	       (1.872046940740317e+04 + s4 *
		(7.446449129977780e+04 + s4 *
		 (2.036310429574584e+05 + s4 *
		  (3.290762690561379e+05 + s4 *
		   (-3.643169244363664e+05 + s4 *
		    (-5.752005633858765e+06 + s4 *
		     (-2.943736461126254e+07 + s4 *
		      (-1.081239626608858e+08 + s4 *
		       (-3.054119954034933e+08 + s4 *
			(-5.778623048717837e+08 + s4 *
			 (1.039349175115788e+08 + s4 *
			  (7.677528988310307e+09 + s4 *
			   (4.726923328634641e+10 + s4 *
			    (2.069296456479017e+11 + s4
			     )))))))))))))))))))))))));
  x12a = 
    - (0.5 / s / s1
       + 0.225 * ls2
       + (3.0 / 112.0) * s1 * ls2
       + 12.0 / 112.0 / s)
    + (1.142857142857143e-01 + s4 *
       (5.392857142857143e-01 + s4 *
	(-1.442321428571429e+00 + s4 *
	 (-1.811702806122448e+00 + s4 *
	  (4.067142325680273e+00 + s4 *
	   (1.814984062161797e+01 + s4 *
	    (1.517123679543115e+01 + s4 *
	     (-9.315442856337989e+01 + s4 *
	      (-4.509730859708386e+02 + s4 *
	       (-1.013717561768030e+03 + s4 *
		(-5.046648273141244e+02 + s4 *
		 (6.729598718538082e+03 + s4 *
		  (3.755561861642751e+04 + s4 *
		   (1.255728420634840e+05 + s4 *
		    (2.790170687704862e+05 + s4 *
		     (1.783022452325151e+05 + s4 *
		      (-2.116256323231557e+06 + s4 *
		       (-1.437886077349634e+07 + s4 *
			(-6.067945424181265e+07 + s4 *
			 (-1.961945489797276e+08 + s4 *
			  (-4.766047893817811e+08 + s4 *
			   (-5.973100212287734e+08 + s4 *
			    (1.987223471390268e+09 + s4 *
			     (1.957094902568444e+10 + s4 *
			      (9.916629255946661e+10
			       ))))))))))))))))))))))))) / s;
  y11a = 
    - ls1 / 6.0
    + 1.0
    + s4 *
    (-1.041666666666666e-01 + s4 *
     (4.830729166666667e-01 + s4 *
      (4.429796006944464e-02 + s4 *
       (-1.968897501627604e+00 + s4 *
	(-2.014702924092610e+00 + s4 *
	 (7.494522048367401e+00 + s4 *
	  (2.972570308847799e+01 + s4 *
	   (4.402886766913980e+01 + s4 *
	    (-4.476216812665643e+01 + s4 *
	     (-5.635225332717891e+02 + s4 *
	      (-2.357775552631778e+03 + s4 *
	       (-6.405099938013067e+03 + s4 *
		(-8.430160130264936e+03 + s4 *
		 (2.624984448061837e+04 + s4 *
		  (2.443423629098088e+05 + s4 *
		   (1.092627425171718e+06 + s4 *
		    (3.496730661033154e+06 + s4 *
		     (7.785502318307757e+06 + s4 *
		      (5.234771883484364e+06 + s4 *
		       (-6.056921196040344e+07 + s4 *
			(-4.314272359664154e+08 + s4 *
			 (-1.932720650644135e+09 + s4 *
			  (-6.836031616103577e+09 + s4 *
			   (-1.957431871170093e+10 + s4 *
			    (-4.100639948144824e+10
			     )))))))))))))))))))))))));
  y12a = 
    - ls2 / 6.0
    + (-8.333333333333337e-02 + s4 *
       (-3.298611111111116e-02 + s4 *
	(-7.272135416666670e-02 + s4 *
	 (6.238635835193449e-01 + s4 *
	  (1.982238628246165e+00 + s4 *
	   (-1.564461411851831e+00 + s4 *
	    (-1.845654074542034e+01 + s4 *
	     (-4.549420266116260e+01 + s4 *
	      (-3.061591747969214e+01 + s4 *
	       (2.356294135776698e+02 + s4 *
		(1.378828319160533e+03 + s4 *
		 (4.589159554966434e+03 + s4 *
		  (9.500180581923400e+03 + s4 *
		   (1.239731247230200e+03 + s4 *
		    (-9.787136158324964e+04 + s4 *
		     (-5.602529184208997e+05 + s4 *
		      (-2.081347283242524e+06 + s4 *
		       (-5.644419377065599e+06 + s4 *
			(-9.152294755832911e+06 + s4 *
			 (1.118954822345829e+07 + s4 *
			  (1.749356103797569e+08 + s4 *
			   (9.354457746353607e+08 + s4 *
			    (3.689108631134369e+09 + s4 *
			     (1.178498366345715e+10 + s4 *
			      (2.967578150982861e+10
			       ))))))))))))))))))))))))) / s;
  y11b = 
    2.0 / 3.0 *
    (-0.25 * ls2
     - 0.125 * s1 * ls2
     - 0.5 / s)
    + (1.333333333333333e+00 + s4 *
       (-7.500000000000000e-01 + s4 *
	(3.472222222223505e-04 + s4 *
	 (-8.611142113095240e-01 + s4 *
	  (-4.453734382750497e+00 + s4 *
	   (-5.206488240046131e+00 + s4 *
	    (1.179395372495230e+01 + s4 *
	     (5.878893157327785e+01 + s4 *
	      (1.076297029275532e+02 + s4 *
	       (-3.409193146082528e+01 + s4 *
		(-1.160440084088181e+03 + s4 *
		 (-5.761051124238334e+03 + s4 *
		  (-1.906396787928903e+04 + s4 *
		   (-4.325100851987218e+04 + s4 *
		    (-3.641117899898082e+04 + s4 *
		     (2.499611696234173e+05 + s4 *
		      (1.736882357980497e+06 + s4 *
		       (6.879687074978024e+06 + s4 *
			(1.938307648180932e+07 + s4 *
			 (3.234298607971404e+07 + s4 *
			  (-4.118373015361849e+07 + s4 *
			   (-6.585704567245671e+08 + s4 *
			    (-3.642893474533921e+09 + s4 *
			     (-1.502089145280693e+10 + s4 *
			      (-5.126202418956630e+10
			       ))))))))))))))))))))))))) / s;
  y12b = 
    - 2.0 / 3.0 *
    (0.25 * ls1
     + 0.125 * s1 * ls1)
    + s4 *
    (0.000000000000000e+00 + s4 *
     (-1.041666666666667e-01 + s4 *
      (-4.752604166666665e-01 + s4 *
       (1.609402126736111e+00 + s4 *
	(6.180223592122394e+00 + s4 *
	 (4.608370399475104e+00 + s4 *
	  (-1.864521149321206e+01 + s4 *
	   (-6.617293552440540e+01 + s4 *
	    (-6.750041678441767e+01 + s4 *
	     (3.421091838234507e+02 + s4 *
	      (2.502962672138401e+03 + s4 *
	       (9.899553935392149e+03 + s4 *
		(2.696438569271003e+04 + s4 *
		 (4.073420721374186e+04 + s4 *
		  (-6.538714673100505e+04 + s4 *
		   (-7.934706454494642e+05 + s4 *
		    (-3.683309834722957e+06 + s4 *
		     (-1.190624155922549e+07 + s4 *
		      (-2.616883483199247e+07 + s4 *
		       (-1.254877523749109e+07 + s4 *
			(2.482848281127249e+08 + s4 *
			 (1.722290365904303e+09 + s4 *
			  (7.869875263182927e+09 + s4 *
			   (2.908162042138783e+10 + s4 *
			    (9.030106198199509e+10
			     )))))))))))))))))))))))));
  x11c = 
    4.0 / 3.0 *
    (0.25 * ls1
     + 0.5 / s * ls2
     + 1.0)
    + s4 *
    (-1.333333333333333e+00 + s4 *
     (-8.888888888888888e-01 + s4 *
      (-8.888888888888886e-02 + s4 *
       (9.523809523809526e-01 + s4 *
	(4.148148148148145e-01 + s4 *
	 (-6.020202020202021e+00 + s4 *
	  (-2.801465201465202e+01 + s4 *
	   (-8.604444444444444e+01 + s4 *
	    (-2.204531590413943e+02 + s4 *
	     (-5.036070175438596e+02 + s4 *
	      (-1.052386724386724e+03 + s4 *
	       (-2.054338164251208e+03 + s4 *
		(-3.981604102564103e+03 + s4 *
		 (-9.099569664902985e+03 + s4 *
		  (-2.980933639846742e+04 + s4 *
		   (-1.282467526881720e+05 + s4 *
		    (-5.865882448009501e+05 + s4 *
		     (-2.591454167195767e+06 + s4 *
		      (-1.085288800189664e+07 + s4 *
		       (-4.326744000683761e+07 + s4 *
			(-1.656019497158343e+08 + s4 *
			 (-6.134249678195915e+08 + s4 *
			  (-2.214286331112399e+09 + s4 *
			   (-7.834394267167847e+09 + s4 *
			    (-2.730646499782967e+10
			     )))))))))))))))))))))))));
  x12c = 
    4.0 / 3.0 *
    (0.25 * ls2
     + 0.5 * ls1 / s
     - 1.0 / s)
    + s4 *
    (-4.444444444444445e-01 + s4 *
     (1.066666666666667e+00 + s4 *
      (2.031746031746032e+00 + s4 *
       (3.407407407407407e+00 + s4 *
	(4.412121212121211e+00 + s4 *
	 (3.008547008547007e+00 + s4 *
	  (-3.974603174603165e+00 + s4 *
	   (-1.474509803921569e+01 + s4 *
	    (1.000389863547753e+01 + s4 *
	     (2.768126984126987e+02 + s4 *
	      (1.612184453227931e+03 + s4 *
	       (6.922702222222224e+03 + s4 *
		(2.566089648622982e+04 + s4 *
		 (8.698871592775041e+04 + s4 *
		  (2.776548014336918e+05 + s4 *
		   (8.498171717171718e+05 + s4 *
		    (2.529593147338937e+06 + s4 *
		     (7.422041001001000e+06 + s4 *
		      (2.178710267026540e+07 + s4 *
		       (6.507617650081301e+07 + s4 *
			(2.012348540007381e+08 + s4 *
			 (6.529391008996632e+08 + s4 *
			  (2.233502252087575e+09 + s4 *
			   (8.011012000036275e+09
			    )))))))))))))))))))))))) / s;
  y11c =
    4.0 / 3.0 *
    (- 0.2 * ls1
     - (47.0 / 250.0) * s1 * ls1
     + 1.0 /*YCf [0]*/)
    + s4 *
    (-2.069333333333333e+00 + s4 *
     (8.719999999999999e-01 + s4 *
      (2.880722222222222e+00 + s4 *
       (3.847295138888887e+00 + s4 *
	(-1.324714843750006e-01 + s4 *
	 (-1.366785869683156e+00 + s4 *
	  (2.912432901570148e+01 + s4 *
	   (1.377373560177257e+02 + s4 *
	    (3.683211673877055e+02 + s4 *
	     (6.589208332005879e+02 + s4 *
	      (3.147176571623555e+02 + s4 *
	       (-3.473293062242155e+03 + s4 *
		(-1.643180504875248e+04 + s4 *
		 (-3.217245593226010e+04 + s4 *
		  (6.816971657848917e+04 + s4 *
		   (9.455403655130938e+05 + s4 *
		    (5.252119204144445e+06 + s4 *
		     (2.156161298611233e+07 + s4 *
		      (7.209517502915230e+07 + s4 *
		       (1.968623175231554e+08 + s4 *
			(3.985861234753939e+08 + s4 *
			 (2.854769767232488e+08 + s4 *
			  (-2.405352201178726e+09 + s4 *
			   (-1.653646328876052e+10 + s4 *
			    (-6.891440613246304e+10
			     )))))))))))))))))))))))));
  y12c =
    4.0 / 3.0 *
    (0.05 * ls2
     + (31.0 / 250.0) * s1 * ls2
     + 62.0 / 125.0 / s)
    +
    (-1.589333333333333e+00 + s4 *
     (2.074666666666666e+00 + s4 *
      (1.307511111111111e+00 + s4 *
       (1.902370238095238e+00 + s4 *
	(-2.753183283730157e-01 + s4 *
	 (-2.011795837607854e+00 + s4 *
	  (2.317028745844202e+01 + s4 *
	   (1.395279499237549e+02 + s4 *
	    (4.617415542216800e+02 + s4 *
	     (1.169356664712524e+03 + s4 *
	      (2.407773115940960e+03 + s4 *
	       (3.787806439753997e+03 + s4 *
		(4.328153450908105e+03 + s4 *
		 (1.237230859426213e+04 + s4 *
		  (1.140755302278274e+05 + s4 *
		   (7.921173630140705e+05 + s4 *
		    (4.140560313559534e+06 + s4 *
		     (1.788189678733328e+07 + s4 *
		      (6.711395833886759e+07 + s4 *
		       (2.248455262052208e+08 + s4 *
			(6.829695936769902e+08 + s4 *
			 (1.906745973315716e+09 + s4 *
			  (5.018118831118233e+09 + s4 *
			   (1.325947551557698e+10 + s4 *
			    (3.974164210617063e+10
			     ))))))))))))))))))))))))) / s;
  x11g = 
    0.5/s1/s
    + (0.225
       +1.392857142857142905e-01*s5)*ls2
    -1.392857142857142905e-01*s
    +
    (-1.028571428571428e+00 + s4 *
     (8.471428571428571e-01 + s4 *
      (2.066887755102041e+00 + s4 *
       (-7.737223639455782e+00 + s4 *
	(-1.894127096861472e+01 + s4 *
	 (2.614173853490259e+01 + s4 *
	  (1.945725675008585e+02 + s4 *
	   (3.577223455759256e+02 + s4 *
	    (-3.000809356095572e+02 + s4 *
	     (-4.032615524189556e+03 + s4 *
	      (-1.452373943396063e+04 + s4 *
	       (-3.072046392481669e+04 + s4 *
		(-1.391312266620251e+04 + s4 *
		 (2.328070089368676e+05 + s4 *
		  (1.398796774203775e+06 + s4 *
		   (5.271362213408929e+06 + s4 *
		    (1.453352873037549e+07 + s4 *
		     (2.467022711079157e+07 + s4 *
		      (-2.235573964279457e+07 + s4 *
		       (-4.256230063687916e+08 + s4 *
			(-2.333919269854781e+09 + s4 *
			 (-9.304322051070610e+09 + s4 *
			  (-2.980794984901221e+10 + s4 *
			   (-7.435449953586630e+10 + s4 *
			    (-1.069561700124279e+11
			     ))))))))))))))))))))))))) / s;
  x12g = 
    -1.0/s1/s3
    + (0.225
       +1.392857142857142905e-01*s5)*ls1
    +1.392857142857142905e-01
    + s4 *
    (-8.785714285714286e-01 + s4 *
     (3.803571428571429e+00 + s4 *
      (-4.224107142857143e+00 + s4 *
       (-1.282816964285714e+01 + s4 *
	(-1.410613839285711e+00 + s4 *
	 (8.748849449936225e+01 + s4 *
	  (2.713393430125957e+02 + s4 *
	   (2.305344956170946e+02 + s4 *
	    (-1.361499873119535e+03 + s4 *
	     (-7.674093481537458e+03 + s4 *
	      (-2.235906157249475e+04 + s4 *
	       (-3.477574377684933e+04 + s4 *
		(4.694713831685032e+04 + s4 *
		 (6.144557079826149e+05 + s4 *
		  (2.876193678217343e+06 + s4 *
		   (9.428248086751513e+06 + s4 *
		    (2.181218902498707e+07 + s4 *
		     (1.960014373454566e+07 + s4 *
		      (-1.370979650344465e+08 + s4 *
		       (-1.057620405406006e+09 + s4 *
			(-4.830789642138412e+09 + s4 *
			 (-1.716500340801855e+10 + s4 *
			  (-4.878161468501544e+10 + s4 *
			   (-9.854983749044119e+10 + s4 *
			    (-2.903370662441536e+10
			     )))))))))))))))))))))))));
  y11g = 
    (one12
     +one24*s5)*ls2
    -one24*s
    +
    (-2.222222222222222e-01 + s4 *
     (-3.555555555555555e-01 + s4 *
      (1.335714285714286e+00 + s4 *
       (5.663029100529116e-02 + s4 *
	(-7.813864557028619e+00 + s4 *
	 (-9.592514358200393e+00 + s4 *
	  (2.885452067709376e+01 + s4 *
	   (1.382276854511959e+02 + s4 *
	    (2.717690198182968e+02 + s4 *
	     (5.523433742799109e+01 + s4 *
	      (-2.117449245378311e+03 + s4 *
	       (-1.134237454353018e+04 + s4 *
		(-3.804944028184589e+04 + s4 *
		 (-8.302121588389931e+04 + s4 *
		  (-3.428715765637223e+04 + s4 *
		   (7.525048372302778e+05 + s4 *
		    (4.729261119399047e+06 + s4 *
		     (1.889666052107433e+07 + s4 *
		      (5.654037027365173e+07 + s4 *
		       (1.161301103584850e+08 + s4 *
			(3.291773313300939e+07 + s4 *
			 (-1.186891924274502e+09 + s4 *
			  (-7.757624122586328e+09 + s4 *
			   (-3.423875879531926e+10 + s4 *
			    (-1.221686895434343e+11
			     ))))))))))))))))))))))))) / s;
  y12g = 
    (one12
     +one24*s5)*ls1
    +one24
    + s4 *
    (2.500000000000000e-01 + s4 *
     (-7.777777777777778e-01 + s4 *
      (-1.319444444444445e-01 + s4 *
       (3.100781250000000e+00 + s4 *
	(8.995925564236112e+00 + s4 *
	 (-1.401699732220365e+00 + s4 *
	  (-6.733832763490223e+01 + s4 *
	   (-2.179581714228348e+02 + s4 *
	    (-3.392540971750742e+02 + s4 *
	     (3.874599704135788e+02 + s4 *
	      (5.242816551916470e+03 + s4 *
	       (2.305925522395688e+04 + s4 *
		(6.566556498627734e+04 + s4 *
		 (1.018677793966759e+05 + s4 *
		  (-1.681607308313629e+05 + s4 *
		   (-2.089027289603405e+06 + s4 *
		    (-1.004095247494084e+07 + s4 *
		     (-3.443335138597512e+07 + s4 *
		      (-8.698335875614998e+07 + s4 *
		       (-1.181238156022712e+08 + s4 *
			(3.097943513985258e+08 + s4 *
			 (3.308313433384167e+09 + s4 *
			  (1.684103848393915e+10 + s4 *
			   (6.582249615623303e+10 + s4 *
			    (2.128777209057203e+11
			     )))))))))))))))))))))))));
  y11h = 
    - (one30
       +g3yh*s5)*ls1
    -g3yh
    + s4 *
    (4.933333333333334e-02 + s4 *
     (-2.311111111111112e-02 + s4 *
      (-2.224000000000000e+00 + s4 *
       (7.232333333333332e-01 + s4 *
	(1.578172986111111e+01 + s4 *
	 (2.416211709449405e+01 + s4 *
	  (-3.959598063732329e+01 + s4 *
	   (-2.759947788659555e+02 + s4 *
	    (-7.380886636688762e+02 + s4 *
	     (-9.932918942272195e+02 + s4 *
	      (1.856073641418639e+03 + s4 *
	       (1.923125423525319e+04 + s4 *
		(8.173048129805524e+04 + s4 *
		 (2.319527552407343e+05 + s4 *
		  (3.717349318895601e+05 + s4 *
		   (-5.222943636352951e+05 + s4 *
		    (-7.233711292342820e+06 + s4 *
		     (-3.611833457725316e+07 + s4 *
		      (-1.293363297186427e+08 + s4 *
		       (-3.510308689950631e+08 + s4 *
			(-6.002790030336890e+08 + s4 *
			 (4.923319619222050e+08 + s4 *
			  (1.021996287753353e+10 + s4 *
			   (5.795577659662904e+10 + s4 *
			    (2.420311559228540e+11
			     )))))))))))))))))))))))));
  y12h = 
    (two15
     +g5yh*s5)*ls2
    -g5yh*s
    +
    (-3.324444444444444e-01 + s4 *
     (1.116266666666667e+00 + s4 *
      (-1.431161904761905e+00 + s4 *
       (-2.013957671957672e+00 + s4 *
	(4.820332070707069e+00 + s4 *
	 (1.856652828161421e+01 + s4 *
	  (6.669558284087865e+00 + s4 *
	   (-1.065326253605862e+02 + s4 *
	    (-4.109576193722526e+02 + s4 *
	     (-8.310775193426895e+02 + s4 *
	      (-1.463401235532169e+02 + s4 *
	       (7.710185473459562e+03 + s4 *
		(4.250505151012226e+04 + s4 *
		 (1.448836751293515e+05 + s4 *
		  (3.212731576151169e+05 + s4 *
		   (1.599895478458714e+05 + s4 *
		    (-2.757129639770347e+06 + s4 *
		     (-1.782153700022562e+07 + s4 *
		      (-7.281778926093897e+07 + s4 *
		       (-2.257081834426685e+08 + s4 *
			(-5.051480832332337e+08 + s4 *
			 (-4.116712930403453e+08 + s4 *
			  (3.468072136925576e+09 + s4 *
			   (2.648970462062139e+10 + s4 *
			    (1.246212962763133e+11
			     ))))))))))))))))))))))))) / s;
  x11m = 
    two3/s1/s3
    - (0.15
       +1.400793650793651091e-01*s5)*ls1
    +9.710317460317460236e-01
    + s4 *
    (-9.865079365079366e-01 + s4 *
     (6.735449735449730e-01 + s4 *
      (1.069986772486772e+01 + s4 *
       (-8.088234126984133e+00 + s4 *
	(-9.488393683862436e+01 + s4 *
	 (-1.104613659178950e+02 + s4 *
	  (4.520050903043505e+02 + s4 *
	   (2.263032330462203e+03 + s4 *
	    (4.528144356865987e+03 + s4 *
	     (-4.977621385472781e+02 + s4 *
	      (-4.134405523596569e+04 + s4 *
	       (-1.912774740214068e+05 + s4 *
		(-5.645364932224272e+05 + s4 *
		 (-1.042899372418414e+06 + s4 *
		  (2.668334026481701e+05 + s4 *
		   (1.302425058380290e+07 + s4 *
		    (7.269199506122653e+07 + s4 *
		     (2.798864998574021e+08 + s4 *
		      (8.294664242644068e+08 + s4 *
		       (1.726827124928086e+09 + s4 *
			(6.854379672929261e+08 + s4 *
			 (-1.663301835146055e+10 + s4 *
			  (-1.128916941039529e+11 + s4 *
			   (-5.115495927470804e+11 + s4 *
			    (-1.875486629722044e+12
			     )))))))))))))))))))))))));
  x12m = 
    one3/s1/s
    + (0.15
       +1.956349206349206338e-01*s5)*ls2
    -1.956349206349206338e-01*s
    +
    (-4.116402116402118e-01 + s4 *
     (3.839576719576720e+00 + s4 *
      (-1.362120181405896e+01 + s4 *
       (6.333374275636177e+00 + s4 *
	(4.400992113596280e+01 + s4 *
	 (8.799205906940279e+01 + s4 *
	  (-7.355329706896123e+01 + s4 *
	   (-9.798142420531591e+02 + s4 *
	    (-3.087794405835338e+03 + s4 *
	     (-4.139714761022082e+03 + s4 *
	      (9.893903743404746e+03 + s4 *
	       (8.727724670576537e+04 + s4 *
		(3.414997993141904e+05 + s4 *
		 (8.873525409640922e+05 + s4 *
		  (1.145022685911650e+06 + s4 *
		   (-3.442151252741758e+06 + s4 *
		    (-3.270692000754736e+07 + s4 *
		     (-1.512427262937129e+08 + s4 *
		      (-5.154009359990964e+08 + s4 *
		       (-1.317309368591246e+09 + s4 *
			(-1.895369621230615e+09 + s4 *
			 (4.105460025633299e+09 + s4 *
			  (4.840340789537030e+10 + s4 *
			   (2.541894467425741e+11 + s4 *
			    (1.019480957855432e+12
			     ))))))))))))))))))))))))) / s;
  y11m = 
    - (two15
       +2.533333333333333270e-02*s5)*ls1
    +1.085777777777777731e+00
    + s4 *
    (-4.826666666666667e-01 + s4 *
     (-9.991111111111112e-01 + s4 *
      (8.401777777777777e+00 + s4 *
       (9.020444444444460e-01 + s4 *
	(-6.298362222222223e+01 + s4 *
	 (-1.078007509920635e+02 + s4 *
	  (1.748204823288690e+02 + s4 *
	   (1.272527054461444e+03 + s4 *
	    (3.383773922149320e+03 + s4 *
	     (4.134469878787685e+03 + s4 *
	      (-1.093786137223340e+04 + s4 *
	       (-9.681588900155625e+04 + s4 *
		(-4.003051832452340e+05 + s4 *
		 (-1.136175884044256e+06 + s4 *
		  (-1.896706134488004e+06 + s4 *
		   (2.055028163695717e+06 + s4 *
		    (3.372409697675041e+07 + s4 *
		     (1.749139685215164e+08 + s4 *
		      (6.494633201615030e+08 + s4 *
		       (1.859672145362082e+09 + s4 *
			(3.641992185669537e+09 + s4 *
			 (2.394752425827114e+08 + s4 *
			  (-4.231528708200645e+10 + s4 *
			   (-2.690766928915524e+11 + s4 *
			    (-1.192878468859270e+12
			     )))))))))))))))))))))))));
  y12m = 
    (one30
     +1.479999999999999927e-01*s5)*ls2
    -1.479999999999999927e-01*s
    +
    (2.613333333333334e-01 + s4 *
     (-2.639822222222223e+00 + s4 *
      (9.003479365079365e+00 + s4 *
       (-1.625396825396826e-02 + s4 *
	(-2.750871380471381e+01 + s4 *
	 (-7.345104875679876e+01 + s4 *
	  (-1.051453357371794e+01 + s4 *
	   (4.870907492174735e+02 + s4 *
	    (1.948371997172542e+03 + s4 *
	     (4.484620023591169e+03 + s4 *
	      (3.208042836140465e+03 + s4 *
	       (-3.051113811024446e+04 + s4 *
		(-1.974025190975534e+05 + s4 *
		 (-7.260294329836837e+05 + s4 *
		  (-1.761322679389347e+06 + s4 *
		   (-1.586166383016621e+06 + s4 *
		    (1.121385854311815e+07 + s4 *
		     (8.353496785294113e+07 + s4 *
		      (3.618230317437328e+08 + s4 *
		       (1.182737707203069e+09 + s4 *
			(2.890655228129366e+09 + s4 *
			 (3.685602051830912e+09 + s4 *
			  (-1.152502310987768e+10 + s4 *
			   (-1.154301052975062e+11 + s4 *
			    (-5.857284377817185e+11
			     ))))))))))))))))))))))))) / s;
  z11m = 
    one12*s5*ls1
    +1.194444444444444420e+00
    + s4 *
    (-1.666666666666667e-01 + s4 *
     (-2.222222222222223e-01 + s4 *
      (-4.444444444444445e-01 + s4 *
       (2.405555555555556e+00 + s4 *
	(4.888888888888884e-01 + s4 *
	 (-1.396031746031746e+01 + s4 *
	  (-3.038095238095239e+01 + s4 *
	   (1.485995370370369e+01 + s4 *
	    (3.122133101851852e+02 + s4 *
	     (1.280971464646465e+03 + s4 *
	      (3.350171085858586e+03 + s4 *
	       (4.447650941506410e+03 + s4 *
		(-1.220601767351572e+04 + s4 *
		 (-1.202936876469494e+05 + s4 *
		  (-5.567766454141866e+05 + s4 *
		   (-1.864567362740959e+06 + s4 *
		    (-4.523666943355881e+06 + s4 *
		     (-5.115481057567388e+06 + s4 *
		      (2.262275962387467e+07 + s4 *
		       (1.992722434468994e+08 + s4 *
			(9.608839982682569e+08 + s4 *
			 (3.600989243544849e+09 + s4 *
			  (1.106733557046591e+10 + s4 *
			   (2.652958447809672e+10 + s4 *
			    (3.440812805312791e+10
			     )))))))))))))))))))))))));
  z12m = 
    one12*s5*ls2
    -one12*s
    +
    (2.222222222222223e-01 + s4 *
     (1.777777777777778e-01 + s4 *
      (-1.917460317460318e+00 + s4 *
       (6.772486772486773e-01 + s4 *
	(1.723905723905724e+00 + s4 *
	 (4.773892773892775e+00 + s4 *
	  (2.007980769230770e+01 + s4 *
	   (5.193120915032680e+01 + s4 *
	    (1.304299965600278e+01 + s4 *
	     (-5.231094402673349e+02 + s4 *
	      (-2.669614089458247e+03 + s4 *
	       (-7.539883371829709e+03 + s4 *
		(-9.456860148258374e+03 + s4 *
		 (3.230499318574015e+04 + s4 *
		  (2.751861525157783e+05 + s4 *
		   (1.148170115793047e+06 + s4 *
		    (3.360716258745590e+06 + s4 *
		     (6.301367820611982e+06 + s4 *
		      (-1.390872809498948e+06 + s4 *
		       (-8.013301441302040e+07 + s4 *
			(-4.650701539897203e+08 + s4 *
			 (-1.901074133659964e+09 + s4 *
			  (-6.274909562234401e+09 + s4 *
			   (-1.666503922249157e+10 + s4 *
			    (-3.008135575916472e+10
			     ))))))))))))))))))))))))) / s;
  
  res [ 0] = x11a;
  res [ 1] = x12a;
  res [ 2] = y11a;
  res [ 3] = y12a;
  res [ 4] = y11b;
  res [ 5] = y12b;
  res [ 6] = x11c;
  res [ 7] = x12c;
  res [ 8] = y11c;
  res [ 9] = y12c;
  res [10] = x11g;
  res [11] = x12g;
  res [12] = y11g;
  res [13] = y12g;
  res [14] = y11h;
  res [15] = y12h;
  res [16] = x11m;
  res [17] = x12m;
  res [18] = y11m;
  res [19] = y12m;
  res [20] = z11m;
  res [21] = z12m;
}