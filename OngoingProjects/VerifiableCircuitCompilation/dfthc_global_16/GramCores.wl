BeginPackage["DFTHCGram`"];
AugmentedCore::usage="AugmentedCore[w,v] returns [w;v^T][w;v^T]^T for real coefficient rows and identity offsets.";
ProjectorGram::usage="ProjectorGram[u,v] returns (u.Transpose[v])^2 elementwise for real row vectors.";
ProjectorGradient::usage="ProjectorGradient[u,w,u0,w0] gives the Euclidean gradient of the squared factor-matrix residual, holding weights fixed.";
SchurRescaleChecks::usage="SchurRescaleChecks[] proves offset-preserving positive core rescaling and convexity of the residual fit.";
OffsetCompressionChecks::usage="OffsetCompressionChecks[] proves offset-preserving core rank reduction and the coupled residual objective.";
GramSymbolicChecks::usage="GramSymbolicChecks[] derives core, Schur, projection, gradient, rank and failure identities exactly.";
Begin["`Private`"];
realQ[x_]:=NumericQ[x] && TrueQ[Im[x]==0] && FreeQ[x,Indeterminate|_DirectedInfinity];
realMatrixQ[x_]:=MatrixQ[x,realQ] && Min[Dimensions[x]]>0;
AugmentedCore[w_?realMatrixQ,v_List] /; VectorQ[v,realQ] && Length[v]==Last[Dimensions[w]]:=
  With[{a=Join[w,{v}]},a.Transpose[a]];
AugmentedCore[___]:=Failure["CoreDomain",<||>];
ProjectorGram[u_?realMatrixQ,v_?realMatrixQ] /; Last[Dimensions[u]]==Last[Dimensions[v]]:=
  (u.Transpose[v])^2;
ProjectorGram[___]:=Failure["ProjectorDomain",<||>];
ProjectorGradient[u_?realMatrixQ,w_?realMatrixQ,u0_?realMatrixQ,w0_?realMatrixQ] /;
  Length[u]==Length[w] && Length[u0]==Length[w0] && Last[Dimensions[u]]==Last[Dimensions[u0]] &&
  Last[Dimensions[w]]==Last[Dimensions[w0]]:=
  4 ((w.Transpose[w] (u.Transpose[u])).u-(w.Transpose[w0] (u.Transpose[u0])).u0);
ProjectorGradient[___]:=Failure["GradientDomain",<||>];

GramSymbolicChecks[]:=Module[
 {a,b,c,d,e,f,g,h,t,x,y,z,q1,q2,tau,u,u0,w,w0,v,aa,k,pp,xx,offset,
  rotation,rotated,singular,principal,schur,schurResult,loss,variables,gradient,expected,
  gram,projection,residual,rankFrame,leftInverse,diagonal,fullRankMatrix,toy,lambda2,
  uh,uth,wh,wth},
 w={{a,b},{c,d}};v={e,f};aa=Join[w,{v}];k=aa.Transpose[aa];
 pp={{1,0},{0,1},{g,h}};xx=pp.w;offset=w.v;
 rotation={{1-t^2,-2t},{2t,1-t^2}}/(1+t^2);rotated=aa.rotation;
 singular={{1,0,q1},{0,0,q2},{q1,q2,tau}};
 principal=Flatten[Table[Det[singular[[#, #]]]& /@ Subsets[Range[3],{j}],{j,3}]];
 schurResult=Reduce[And@@Thread[principal>=0],{q1,q2,tau},Reals];
 schur=Resolve[ForAll[{q1,q2,tau},Equivalent[And@@Thread[principal>=0],q2==0 && tau>=q1^2]],Reals];
 u=Array[uh,{2,3}];u0=Array[uth,{3,3}];
 w=Array[wh,{2,2}];w0=Array[wth,{3,2}];
 loss=Total[Flatten[w.Transpose[w] (u.Transpose[u])^2]]-
   2 Total[Flatten[w.Transpose[w0] (u.Transpose[u0])^2]];
 variables=Flatten[u];gradient=ArrayReshape[D[loss,{variables}],Dimensions[u]];
 expected=4 ((w.Transpose[w] (u.Transpose[u])).u-(w.Transpose[w0] (u.Transpose[u0])).u0);
 gram={{a,b},{b,c}};projection={1,b/a};
 residual=gram-Outer[Times,gram[[All,1]],gram[[1]]]/a;
 rankFrame={{1,0},{0,1},{a,b}};leftInverse={{1,0,0},{0,1,0}};
 diagonal=DiagonalMatrix[{c,d}];fullRankMatrix=rankFrame.diagonal.Transpose[rankFrame];
 toy={{1,0},{0,1},{-1,-1}};
 lambda2=Function[frame,Total[(Total[Abs[Most[frame]]]+Abs[Last[frame]+Total[Most[frame]]])^2]/4];
 <|"schur_condition"->StringReplace[ToString[schurResult,InputForm],RegularExpression["DFTHCGram`Private`|\\$[0-9]+"]->""],
   "checks"-><|
     "core_blocks"->(Take[k,2,2]===Take[aa,2].Transpose[Take[aa,2]] && k[[1;;2,3]]===offset),
     "two_body_from_core"->TrueQ[And@@Thread[Flatten[Map[Expand,xx.Transpose[xx]-pp.k[[1;;2,1;;2]].Transpose[pp],{2}]]==0]],
     "offset_contraction"->TrueQ[And@@Thread[Expand[xx.v-pp.k[[1;;2,3]]]==0]],
     "offset_norm"->TrueQ[Expand[k[[3,3]]-v.v]==0],
     "orthogonal_gauge"->TrueQ[And@@Thread[Flatten[Factor[rotated.Transpose[rotated]-k]]==0]],
     "singular_schur_range_condition"->TrueQ[schur],
     "schur_infeasible_out_of_range"->TrueQ[Reduce[q2!=0 && And@@Thread[principal>=0],{q1,q2,tau},Reals]===False],
     "core_psd_quadratic_form"->TrueQ[Expand[{x,y,z}.k.{x,y,z}-Total[(Transpose[aa].{x,y,z})^2]]==0],
     "core_rank_exact_case"->(MatrixRank[AugmentedCore[{{1,0},{0,1}},{-1,-1}]]===2),
     "normal_equations"->TrueQ[Expand[gram[[1]].{x-projection[[1]](x+b y/a),y}-0]==0],
     "schur_projection_error"->TrueQ[Expand[{x,y}.residual.{x,y}-y^2(c-b^2/a)]==0],
     "schur_projection_nonnegative"->TrueQ[Resolve[ForAll[{a,b,c,y},a>0 && a c>=b^2,y^2(a c-b^2)>=0],Reals]],
     "projector_gradient"->TrueQ[And@@Thread[Flatten[Map[Expand,gradient-expected,{2}]]==0]],
     "gradient_cross_weights_not_specialized"->AllTrue[Diagonal[w.Transpose[w0]],!TrueQ[#==0]&],
     "left_inverse"->(leftInverse.rankFrame===IdentityMatrix[2]),
     "rank_obstruction_minor"->TrueQ[Expand[Det[fullRankMatrix[[1;;2,1;;2]]]-c d]==0],
     "rank_obstruction_support"->TrueQ[And@@Thread[Flatten[Expand[leftInverse.fullRankMatrix.Transpose[leftInverse]-diagonal]]==0]],
     "normalization_not_core_invariant"->(lambda2[toy]===1/2 && FullSimplify[lambda2[toy.{{1,1},{-1,1}}/Sqrt[2]]]===1),
     "noncommuting_projectors"->With[{p1={{1,0},{0,0}},p2=Outer[Times,{3/5,4/5},{3/5,4/5}]},p1.p2=!=p2.p1],
     "identity_limit"->TrueQ[Limit[rotation,t->0]===IdentityMatrix[2]],
     "vanishing_projection_error"->TrueQ[Limit[y^2(c-b^2/a),c->b^2/a,Assumptions->a>0]==0],
     "invalid_core_rejected"->FailureQ[AugmentedCore[{{1,2}},{1}]],
     "invalid_projectors_rejected"->FailureQ[ProjectorGram[{},{{1}}]],
     "invalid_gradient_rejected"->FailureQ[ProjectorGradient[{{1}},{{1,2}},{{1}},{{1}}]]|>|>];
OffsetCompressionChecks[]:=Module[
 {a,b,c,d,f,g,h,s1,s2,y11,y12,y21,y22,oo,vv,ww,zz,t,
  o,v,w,z,q,ideal,vars,
  reduced,offsetError,normError,projectorError,rankError,e={{a,b},{b,c}},
  x={{d,f},{g,h}},s={s1,s2},y,delta,quartic,derivative,gradient,
  ys={{y11,y12},{y21,y22}},coupled,coupledGradient},
 o=Array[oo,3];v=Array[vv,3];w=Array[ww,{2,3}];z=Array[zz,3];
 q=IdentityMatrix[3]-Outer[Times,o,o];vars=Join[o,v];
 ideal=GroebnerBasis[{o.o-1,o.v},vars];
 reduced=Function[expression,Last[PolynomialReduce[Expand[expression],ideal,vars]]];
 offsetError=reduced /@ (w.q.v-w.v);
 normError=reduced[v.q.v-v.v];
 projectorError=reduced[z.q.z-Total[(q.z)^2]];
 rankError=Last[PolynomialReduce[Expand[CharacteristicPolynomial[q,t]+t(t-1)^2],
   GroebnerBasis[{o.o-1},o],o]];
 y=Transpose[x].s;delta=Tr[(e-Outer[Times,y,y]).(e-Outer[Times,y,y])]-Tr[e.e];
 quartic=(s.x.Transpose[x].s)^2-2 s.x.e.Transpose[x].s;
 derivative=D[quartic,{s}];gradient=4 ((s.x.Transpose[x].s) x.Transpose[x]-x.e.Transpose[x]).s;
 coupled=Total[Flatten[(ys.Transpose[ys])^2]]-2 Tr[ys.e.Transpose[ys]];
 coupledGradient=ArrayReshape[D[coupled,{Flatten[ys]}],{2,2}];
 <|"offset_contraction_preserved"->AllTrue[offsetError,#===0&],
   "offset_norm_preserved"->(normError===0),
   "psd_orthogonal_projector"->(projectorError===0),
   "copy_rank_drops_by_one"->(rankError===0),
   "rank_one_core_removed"->AllTrue[Flatten[Map[Expand,
    w.Transpose[w]-w.q.Transpose[w]-Outer[Times,w.o,w.o],{2}]],#===0&],
   "residual_quartic_identity"->(Expand[delta-quartic]===0),
   "residual_quartic_gradient"->AllTrue[Expand[derivative-gradient],#===0&],
   "coupled_gradient"->AllTrue[Flatten[Map[Expand,coupledGradient-4(ys.Transpose[ys].ys-ys.e),{2}]],#===0&]|>];
SchurRescaleChecks[]:=Module[
 {o=Array[oo,3],v=Array[vv,3],w=Array[ww,{2,3}],q,pv,a,vars,ideal,reduce,
  s,t1,t2,z1,z2,h1,h2,target,loss,gram,hessian},
 vars=Join[o,v];ideal=GroebnerBasis[{o.o-1,o.v},vars];
 reduce=Function[expression,Last[PolynomialReduce[
   Numerator[Together[expression]],ideal,vars]]];
 q=IdentityMatrix[3]-Outer[Times,o,o];pv=Outer[Times,v,v]/(v.v);
 a=pv+s(IdentityMatrix[3]-pv);
 h1={{1,0},{0,0}};h2={{1,1},{1,1}};target={{2,3},{3,5}};
 loss=Total[Flatten[(target-t1 h1-t2 h2)^2]];
 gram=Outer[Tr[#1.#2]&,{h1,h2},{h1,h2},1];
 hessian=D[loss,{{t1,t2},2}];
 <|"schur_rescale_core_identity"->AllTrue[Flatten[Map[reduce,
    w.a.q.Transpose[a].Transpose[w]-w.pv.Transpose[w]-s^2 w.(q-pv).Transpose[w],{2}]],#===0&],
   "schur_rescale_offset_preserved"->AllTrue[reduce/@(w.a.q.v-w.v),#===0&],
   "schur_rescale_norm_preserved"->(reduce[v.q.v-v.v]===0),
   "scale_fit_hessian"->(hessian===2 gram),
   "scale_fit_psd_gram"->TrueQ[Expand[{z1,z2}.gram.{z1,z2}-Total[Flatten[(z1 h1+z2 h2)^2]]]==0],
   "nonnegative_fit_boundary"->(First[Minimize[{(t1+1)^2+(t2-2)^2,t1>=0 && t2>=0},{t1,t2},Reals]]===1)|>];
End[];
EndPackage[];
