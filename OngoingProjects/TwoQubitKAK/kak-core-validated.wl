X=PauliMatrix[1];Y=PauliMatrix[2];Z=PauliMatrix[3];
kp=KroneckerProduct;XX=kp[X,X];YY=kp[Y,Y];ZZ=kp[Z,Z];
MB=(1/Sqrt[2]){{1,0,0,I},{0,I,1,0},{0,I,-1,0},{1,0,0,-I}};MBd=ConjugateTranspose[MB];
S=Transpose[Chop@Re@{Diagonal[MBd.XX.MB],Diagonal[MBd.YY.MB],Diagonal[MBd.ZZ.MB]}];
can[c_]:=MatrixExp[I(c[[1]]XX+c[[2]]YY+c[[3]]ZZ)];
factorKron[L_]:=Module[{R,u,s,v},
  R=ArrayReshape[Table[L[[2(i1-1)+i2,2(j1-1)+j2]],{i1,2},{j1,2},{i2,2},{j2,2}],{4,4}];
  {u,s,v}=SingularValueDecomposition[R];
  {ArrayReshape[Sqrt[s[[1,1]]]u[[All,1]],{2,2}],ArrayReshape[Sqrt[s[[1,1]]]Conjugate[v[[All,1]]],{2,2}]}];

kak[U_]:=Module[{gp,u,uM,M2,t=Sqrt[2.],es,P,theta,K1,K2,c,corr,Lleft,Rright,A,B},
  gp=Det[U]^(1/4);u=U/gp;uM=MBd.u.MB;M2=Transpose[uM].uM;
  es=Eigensystem[Re[M2]+t Im[M2]];P=Transpose[es[[2]]];
  If[Re[Det[P]]<0,P[[All,1]]=-P[[All,1]]];K2=P;
  theta=Arg[Diagonal[Transpose[K2].M2.K2]]/2;
  K1=uM.K2.DiagonalMatrix[Exp[-I theta]];
  If[Re[Det[K1]]<0,K1[[All,1]]=-K1[[All,1]];theta[[1]]+=Pi];
  c=Transpose[S].theta/4;                                  (* orthogonal-column projection *)
  corr=MB.DiagonalMatrix[Exp[I(theta-S.c)]].MBd;           (* local residual (diag +-1 in magic) *)
  Lleft=MB.K1.MBd.corr;  Rright=MB.Transpose[K2].MBd;      (* fold corr into left local *)
  A=factorKron[Lleft]; B=factorKron[Rright];
  <|"gp"->gp,"c"->c,"A"->A,"B"->B,"corrLocalErr"->Chop@Norm[corr-kp@@factorKron[corr],"Frobenius"]|>];

SeedRandom[7];
ab=kp[MatrixExp[I 0.4(0.3X+0.7Y+0.2Z)],MatrixExp[I 0.9(X-0.5Z)]];
CNOT={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};CZ=DiagonalMatrix[{1,1,1,-1}];
SWAP={{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};iSWAP={{1,0,0,0},{0,0,I,0},{0,I,0,0},{0,0,0,1}};
DCX={{1,0,0,0},{0,0,0,1},{0,1,0,0},{0,0,1,0}};cphase=DiagonalMatrix[{1,1,1,Exp[I 0.7]}];
battery=<|"Identity"->IdentityMatrix[4],"AxB"->ab,"CNOT"->CNOT,"CZ"->CZ,"SWAP"->SWAP,
  "iSWAP"->iSWAP,"DCX"->DCX,"cphase"->cphase,
  "r1"->RandomVariate[CircularUnitaryMatrixDistribution[4]],
  "r2"->RandomVariate[CircularUnitaryMatrixDistribution[4]],
  "r3"->RandomVariate[CircularUnitaryMatrixDistribution[4]],
  "r4"->RandomVariate[CircularUnitaryMatrixDistribution[4]]|>;

rep=KeyValueMap[Function[{name,U},
  Module[{k,full},
    k=kak[N[U]];
    full=Chop@Norm[k["gp"](kp@@k["A"]).can[k["c"]].(kp@@k["B"])-N[U],"Frobenius"];
    {name,NumberForm[k["c"]/(Pi/4),3],k["corrLocalErr"],full}]],battery];
Print[Grid[Prepend[rep,{"gate","coords/(Pi/4)","corrLocalErr","fullReconErr"}],Frame->All,Alignment->Left]];
