(* Created with the Wolfram Language : www.wolfram.com *)
<|"S symbolic all pass?" -> True, "A eta=(gci+gba)/(2gd)" -> 1., 
 "A dep-(1/T2+gd)" -> 7.275957614183426*^-12, 
 "A Lindblad DSolve vs Liouvillian-exp |diff|" -> 1.7763568394002505*^-15, 
 "A Lindblad z(tf)" -> 0.6401926858764212, 
 "A Ito meanZ(tf) [nSteps=500,N=1000]" -> 0.6423692600161793, 
 "A |Ito-Lindblad| < 4/sqrt(N)?" -> True, 
 "A det stability threshold nSteps*" -> 4737.410112522892, 
 "A det EM noProj z(tf) @1000 (below thr, blows up)" -> 4.725220317725831, 
 "A det EM proj REL ERROR @{8000,16000,32000} (O(dt), shrinking, NOT \
converged)" -> {0.27681530375559465, 0.13037996404792512, 
   0.06332361767174663}, "A det EM rel error monotonically shrinks?" -> True, 
 "B {Lindblad, EMproj, Ito} meanZ(tf)" -> {0.1764842431860447, 
   0.1812749694033881, 0.1918432308760084}, 
 "B EM & Ito reproduce Lindblad < 4/sqrt(N)?" -> True, 
 "B EMproj maxr <= 1 (r<=1 every path)?" -> True, 
 "B EMnoProj overshoots r>1?" -> True, 
 "D Std[z(tf)] {Ito, EMproj, half-diffusion}" -> 
  {0.6147965680520422, 0.5992966491259936, 0.376145916756089}, 
 "D EM reproduces Ito 2nd moment AND half-diffusion caught (teeth)?" -> True, 
 "E1 equatorial: EM & Ito reproduce Lindblad, EMproj r<=1?" -> True, 
 "E2 {eta, Ito mean final Bloch radius (mixed<1)}" -> 
  {0.3333333333333333, 0.5577266511115374}, 
 "E2 eta<1: EM & Ito reproduce Lindblad; mixed r<1; EMproj r<=1?" -> True, 
 "G discriminant>0 (overdamped/Zeno, no oscillation)?" -> True, 
 "G exact slow eigenvalue vs leading Zeno -(g1+wx^2/dep)" -> 
  {-0.9291451890553151, -0.8902439024390244}, 
 "G strong-meas ensemble reproduces Lindblad, EMproj r<=1?" -> True, 
 "ALL SYMBOLIC PASS" -> True, "ALL NUMERIC PASS" -> True|>
