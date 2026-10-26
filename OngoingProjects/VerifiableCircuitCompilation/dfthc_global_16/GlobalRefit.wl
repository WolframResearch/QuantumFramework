BeginPackage["DFTHCGlobal`"];
GlobalSymbolicChecks::usage="GlobalSymbolicChecks[] derives global Gram gradients, tight-frame offset preservation, polar pullbacks and fixed-frame core updates with symbolic coefficients.";
Begin["`Private`"];
GlobalSymbolicChecks[]:=Module[
 {xh,th,uh,wh,fh,lh,gh,dh,sh,ah,oh,kh,hh,
  x,t,loss,gradient,expected,u,w,p,blocks,allx,residual,gx,gw,gu,
  frame,l,g,d,h,ideal,vars,reduce,projector,a,omega,k,df,dz,
  tangent,s,radial,core,coreLoss,compressed,whitening,alpha,probability,
  gamma,tau,ga,allocationLoss,allocationGradient,allocationExpected,
  chamberWeights,chamberSigns,chamberOffset,chamberIdentitySign,chamberNorm,
  chamberLoss,chamberGradient,chamberExpected,skew1,skew2,symmetricH,
  aa,bb,cc},
 x=Array[xh,{3,2}];t=Array[th,{3,3}];t=(t+Transpose[t])/2;
 residual=x.Transpose[x]-t;loss=Total[Flatten[residual^2]];
 gradient=ArrayReshape[D[loss,{Flatten[x]}],Dimensions[x]];
 expected=4 residual.x;
 u=Array[uh,{2,2,2}];w=Array[wh,{2,2,2}];
 p=Map[Function[rows,Transpose[{rows[[All,1]]^2,
    Sqrt[2] rows[[All,1]] rows[[All,2]],rows[[All,2]]^2}]],u];
 blocks=MapThread[Transpose[#1].#2&,{p,w}];
 allx=Join[blocks[[1]],blocks[[2]],2];
 residual=allx.Transpose[allx]-t;
 gx=Map[Take[4 residual.allx,All,#]&,{{1,2},{3,4}}];
 gw=MapThread[#1.#2&,{p,gx}];
 gu=MapThread[Function[{rows,weights,derivatives},With[
   {packed=weights.Transpose[derivatives]},
   MapThread[Function[{row,v},2 {{v[[1]],v[[2]]/Sqrt[2]},
      {v[[2]]/Sqrt[2],v[[3]]}}.row],{rows,packed}]]],{u,w,gx}];
 loss=Total[Flatten[residual^2]];
 frame=Array[fh,{3,2}];l=Array[lh,{2,2}];
 g=Array[gh,{3,2}];d=Array[dh,{3,2}];vars=Flatten[frame];
 ideal=GroebnerBasis[Flatten[Transpose[frame].frame-IdentityMatrix[2]],vars];
 reduce=Function[expression,Last[PolynomialReduce[
   Numerator[Together[expression]],ideal,vars]]];
 s=frame.l;projector=IdentityMatrix[3]-frame.Transpose[frame];
 tangent=g-frame.(Transpose[frame].g+Transpose[g].frame)/2;
 h=DiagonalMatrix[{hh[1],hh[2]}];a=Transpose[frame].g;
 omega=(Transpose[frame].d-Transpose[d].frame)/
   {{2 hh[1],hh[1]+hh[2]},{hh[1]+hh[2],2 hh[2]}};
 k=(a-Transpose[a])/{{2 hh[1],hh[1]+hh[2]},{hh[1]+hh[2],2 hh[2]}};
 df=projector.d.Inverse[h]+frame.omega;
 dz=projector.g.Inverse[h]+frame.k;
 radial=Array[sh,{2,2}];radial=(radial+Transpose[radial])/2;
 core=Array[ah,{2,2}];core=(core+Transpose[core])/2;
 compressed=Transpose[frame].t.frame;
 coreLoss=Total[Flatten[(frame.core.Transpose[frame]-t)^2]];
 whitening=Total[Flatten[(core-compressed)^2]]+Total[Flatten[t^2]]-
   Total[Flatten[compressed^2]];
 gamma=Array[hh,3];probability=Exp[gamma]/Total[Exp[gamma]];
 alpha=Sqrt[tau probability];ga=Array[gh,3];allocationLoss=ga.alpha;
 allocationGradient=D[allocationLoss,{gamma}];
 allocationExpected=(ga alpha-probability (ga.alpha))/2;
 chamberWeights=Array[wh,3];chamberSigns=Array[sh,3];
 chamberNorm=chamberSigns.chamberWeights+
   chamberIdentitySign(chamberOffset+Total[chamberWeights]);
 chamberLoss=chamberNorm^2/4;
 chamberGradient=D[chamberLoss,{chamberWeights}];
 chamberExpected=chamberNorm(chamberSigns+chamberIdentitySign)/2;
 skew1={{0,oh[1],oh[2]},{-oh[1],0,oh[3]},{-oh[2],-oh[3],0}};
 skew2={{0,kh[1],kh[2]},{-kh[1],0,kh[3]},{-kh[2],-kh[3],0}};
 symmetricH=Array[hh,{3,3}];symmetricH=(symmetricH+Transpose[symmetricH])/2;
 <|"global_gram_gradient"->AllTrue[Flatten[Map[Expand,gradient-expected,{2}]],#===0&],
   "joint_projector_gradient"->AllTrue[Expand[D[loss,{Flatten[u]}]-Flatten[gu]],#===0&],
   "joint_weight_gradient"->AllTrue[Expand[D[loss,{Flatten[w]}]-Flatten[gw]],#===0&],
   "tight_frame_offset_identity"->AllTrue[reduce/@Flatten[Transpose[s].s-Transpose[l].l],#===0&],
   "stiefel_tangent_identity"->AllTrue[reduce/@Flatten[Transpose[frame].tangent+Transpose[tangent].frame],#===0&],
   "polar_pullback_adjoint"->(reduce[Total[Flatten[g df]]-Total[Flatten[dz d]]]===0),
   "polar_radial_limit"->AllTrue[reduce/@Flatten[projector.frame.radial.Inverse[h]],#===0&],
   "sylvester_adjoint_identity"->TrueQ[Expand[
     Tr[Transpose[symmetricH.skew1+skew1.symmetricH].skew2]-
     Tr[Transpose[skew1].(symmetricH.skew2+skew2.symmetricH)]]==0],
   "fixed_frame_core_projection"->(reduce[coreLoss-whitening]===0),
   "positive_allocation_preserves_norm"->TrueQ[FullSimplify[Total[alpha^2]-tau,
     Assumptions->tau>0 && Element[gamma,Reals]]==0],
   "allocation_gradient"->AllTrue[FullSimplify[allocationGradient-allocationExpected,
     Assumptions->tau>0 && Element[gamma,Reals]],#===0&],
   "lcu_sign_chamber_gradient"->AllTrue[Expand[chamberGradient-chamberExpected],#===0&],
   "psd_rank_one_boundary"->(First[Minimize[
     {(aa-3)^2+(bb-2)^2+(cc+1)^2,aa>=0 && bb>=0 && cc>=0 &&
      aa bb==0 && aa cc==0 && bb cc==0},{aa,bb,cc},Reals]]===5),
   "negative_offset_fixture"->PositiveDefiniteMatrixQ[{{2,-1},{-1,3}}],
   "unsupported_positive_offset_fixture"->!PositiveDefiniteMatrixQ[-IdentityMatrix[2]]|>];
End[];
EndPackage[];
