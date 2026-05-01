(* ::Package:: *)

InterThermalAverages[]:=
Module[{dataC0=Import["files/dataC0.dat"],
dataC1=Import["files/dataC1V2.dat"],
dataC2=Import["files/dataC2V2.dat"],
m0,b0,m1,b1,m2,b2,facC0=1/(2(2Pi)^4),fac1=1/(8(2Pi)^4) 2,factor=(81 Sqrt[3])/(16384 \[Pi]^8),
dataC00={},
dataC01={},
dataC02={}
},
Do[AppendTo[dataC00,{El[[1]],El[[2]](1/(2!*3!))facC0(*(El[[1]]*2Pi^2)^3/(El[[1]](BesselK[2,El[[1]]])^3)*)(*mphi(9mphi^2)/mphi^9(0.1/mphi)^23mphi^3*)}],{El,dataC0}];
Do[AppendTo[dataC01,{El[[1]],El[[2]](El[[1]]/3 *1/3!)fac1(*(El[[1]]*2Pi^2)^3/(BesselK[2,El[[1]]])^3*)(*1/mphi(9mphi^2)/mphi^9(0.1/mphi)^2(3mphi^2)(3mphi^3)*)}],{El,dataC1}];
Do[AppendTo[dataC02,{El[[1]],El[[2]](El[[1]]/3*1/2! 1/2!)factor(*(El[[1]]*2Pi^2)^3/( BesselK[2,El[[1]]])^3*)}],{El,dataC2}];
averagedCi=Interpolation[dataC00/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];
averagedC1i=Interpolation[dataC01/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];(*thermal average with p^2/E outside the cross section*)
averagedC2i=Interpolation[dataC02/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];
m0=(averagedCi[dataC0[[1]][[1]]]-averagedCi[dataC0[[2]][[1]]])/(Log10[dataC0[[1]][[1]]]-Log10[dataC0[[2]][[1]]]);
b0=10^(averagedCi[dataC0[[1]][[1]]]-m0 Log10[dataC0[[1]][[1]]]);
m2=(averagedC2i[dataC2[[1]][[1]]]-averagedC2i[dataC2[[2]][[1]]])/(Log10[dataC2[[1]][[1]]]-Log10[dataC2[[2]][[1]]]);
b2=10^(averagedC2i[dataC2[[1]][[1]]]-m2 Log10[dataC2[[1]][[1]]]);
m1=(averagedC1i[dataC1[[1]][[1]]]-averagedC1i[dataC1[[2]][[1]]])/(Log10[dataC1[[1]][[1]]]-Log10[dataC1[[2]][[1]]]);
b1=10^(averagedC1i[dataC1[[1]][[1]]]-m1 Log10[dataC1[[1]][[1]]]);
averagedCif[xphi_?NumericQ,mphi_?NumericQ]:=0.27/mphi^5 Piecewise[{{b0 xphi^(m0-1) (xphi*2Pi^2)^3/( BesselK[2,xphi])^3,xphi<10^-6},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedCi[xphi]/xphi ,xphi>=10^-6}}];
averagedC1if[xphi_?NumericQ,mphi_?NumericQ]:=0.81/mphi^5 Piecewise[{{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 b1 xphi^m1,xphi<10^-2},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedC1i[xphi] ,xphi>=10^-2}}];
averagedC2if[xphi_?NumericQ,mphi_?NumericQ]:=1/mphi^5 Piecewise[{{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 b2 xphi^m2 ,xphi<10^-2},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedC2i[xphi],10^-2<=xphi}}];
C2cannibal[xphi_,mphi_]:=HeavisideTheta[xphi-0.5](averagedC1if[xphi,mphi]-averagedC2if[xphi,mphi])];



datag = Import["files/gEnergy.dat"];
datagentropy = Import["files/gEntropy.dat"];
For[i =0,i<=Length[datag]-1,i++;datag[[i]][[1]]=Log10[datag[[i]][[1]]]];
gg=Interpolation[datag,Method->"Spline"];
grell[T_]:=  gg[Log10[T]];
grel[T_]:=Piecewise[{{grell[T],T<=10^4},{grell[10^4],T>10^4}}]
(*grelx[x_,m_]:=grel[m/x];*)
(*For[i =0,i<=Length[datagentropy]-1,i++;datagentropy[[i]][[1]]=Log10[datagentropy[[i]][[1]]]];*)
gentropy=Interpolation[datagentropy/.{x_,z_}->{Log10[x],z},Method->"Spline"];
gent[xx_,mphi_]:=gentropy[Log10[mphi/xx]];
gentTemp[T_]:=gentropy[Log10[T]];
Clear[s];
s[T_]:=Piecewise[{{(2*Pi^2 gentropy[Log10[T]](*grel[T]*)*T^3)/45,T<1000},{(2*Pi^2 gentropy[Log10[1000]](*grel[1000]*)*T^3)/45,T>=1000}}];
\[Rho]eqq[T_?NumericQ]:=grel[T] Pi^2/30 T^4;
H[T_?NumericQ]:=Module[{Mpl=2.4*10^18},Sqrt[\[Rho]eqq[T]/(3Mpl^2)]];(*(1.66*Sqrt[grel[T]](*10*))/(1.22089*10^19)T^2;*)
Hbar[T_?NumericQ]:=H[T]/(1+T/(3gentTemp[T])*gentTemp'[T]);
ssc[x_,mphi_]:=(*Piecewise[{{*)(2*Pi^2*(*gent[x,mphi]*)106.69594335689291(mphi/x)^3)/45;(*,x>=mphi/1000},{(2*Pi^2*gent[0.1/1000,mphi](mphi/x)^3)/45,x<mphi/1000}}];*)
nphieq=Compile[{{x,_Real,0},{m,_Real,0}},
1/(x*2Pi^2)*m^3*BesselK[2,x],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True];
\[Rho]eq[T_,m_]:=1/(2Pi^2) m^3 T(BesselK[1,m/T]+3 T/m BesselK[2,m/T]);


neq[T_,m_]:=(m^2 T BesselK[2,m/T])/(2 \[Pi]^2)


Yobs[mdm_]:= 9.711808033846503`*^-48 /(mdm(2*43*Pi^2/(11*45) ((8.617333262145`*^-14)2.72548 )^3))


D\[Rho]eq[T_,m_]:=(m ((m^3/T+12 m T) BesselK[0,m/T]+(5 m^2+24 T^2) BesselK[1,m/T]))/(2 \[Pi]^2);


(*to calculate the relic density, remember that 0.12 = \[CapitalOmega]h^2= (mphi s Yobs)/\[Rho]c, Hence Yobs = (9.711808033846503`*^-48 gev^4)/(mphi(2*43*Pi^2/(11*45)((8.617333262145`*^-14)2.72548 gev)^3))*)


Clear[\[Rho]SM];
\[Rho]SM[T_?NumericQ,Trh_?NumericQ]:=grel[Trh] Pi^2/30 T^4;


Peq[T_,m_]:=T (m^2 T BesselK[2,m/T])/(2 \[Pi]^2)


DD\[Rho]eq[T_,m_]:=(m (m^2+9 T^2) (4 m T BesselK[0,m/T]+(m^2+8 T^2) BesselK[1,m/T]))/(2 \[Pi]^2 T^3)


Neq[T_?NumericQ,m_?NumericQ,a_?NumericQ]:=a^3 (m^2 T BesselK[2,m/T])/(2 \[Pi]^2);(*number in equilibrium.*)


Clear[EnergyDensities]
EnergyDensities[]:=Module[{Mpl = 2.4*10^18(*GeV*),ai = 1,
af=10^17,arh,Trh=150,Hsol,TsoldatStandard,TStandardI,afinal,Tnonst,domain},
TsoldatStandard = Table[{ ai (s[Trh]/s[10^T])^(1/3) ,10^T},{T,Log10[10^-13],Log10[1.Trh],(Log10[1.Trh]-Log10[10^-13])/10000}];
TStandardI=Interpolation[TsoldatStandard,InterpolationOrder->2,Method->"Spline"];
domain= TStandardI(*[[1]][[1]][[2]]*)["Domain"](*[[1]][[2]]*);
afinal = af;
Hsol[a_]:=Sqrt[\[Rho]eqq[TStandardI[a]]/(3Mpl^2)](*Total Hubble rate*);
Return[{Hsol,TStandardI,ai,afinal,domain(*,TsoldatStandard*)}]
];


en=EnergyDensities[];


Tstandard=en[[2]];


Clear[\[CapitalGamma]htophiphii]
\[CapitalGamma]htophiphii=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
HeavisideTheta[1-(4mphi^2)/mh^2](v^2 mh)/(16Pi^3) Sqrt[1-(4mphi^2)/mh^2]T BesselK[1,mh/T]
](*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
];
Clear[\[CapitalGamma]htophiphii2]
\[CapitalGamma]htophiphii2=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
HeavisideTheta[1-(4mphi^2)/mh^2] (v^2 mh^2)/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]  (T BesselK[2,mh/T])
](*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
];
\[CapitalGamma]htophiphi[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphii[T,mphi];
\[CapitalGamma]htophiphi2[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphii2[T,mphi];


\[CapitalGamma]htophiphi2Energydens=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,(*mh=125*)mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
HeavisideTheta[1-(4mphi^2)/mh^2] (v^2 mh^2)/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]  (T BesselK[2,mh/T])
]
,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


\[CapitalGamma]htophiphii2Energydens[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphi2Energydens[T,mphi];


Clear[\[CapitalGamma]\[Phi]to2chi]
\[CapitalGamma]\[Phi]to2chi=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(*Piecewise[{{*)coupl mphi/(16Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]T BesselK[1,mphi/T](*,T<=Tc},{0,T>Tc}}]*)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


Clear[\[CapitalGamma]\[Phi]to2chi2];
\[CapitalGamma]\[Phi]to2chi2=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(*Piecewise[{{*)coupl 1/(32Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]mphi^2 (T BesselK[2,mphi/T]-E^(-(mphi/T)) Sqrt[1/mphi] Sqrt[\[Pi]/2] T^(3/2)(*-mphi 10^p2[[2]][mphi/T]*))(*,T<=Tc},{0,T>Tc}}]*)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


Clear[\[CapitalGamma]\[Phi]to2chiEnergy];
\[CapitalGamma]\[Phi]to2chiEnergy=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
coupl 1/(32Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]mphi^2 (T BesselK[2,mphi/T])
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


\[CapitalGamma]\[Phi]to2\[Chi][T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chi[T,mphi,mchi];


\[CapitalGamma]\[Phi]to2\[Chi]C2[T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chi2[T,mphi,mchi];


\[CapitalGamma]\[Phi]to2chiEner[T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chiEnergy[T,mphi,mchi];


Clear[D\[CapitalGamma]\[Phi]to2chiEnergy]
D\[CapitalGamma]\[Phi]to2chiEnergy=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(coupl Sqrt[1-(4 mchi^2)/mphi^2] mphi^2 (mphi BesselK[1,mphi/T]+3 T BesselK[2,mphi/T]) HeavisideTheta[1-(4 mchi^2)/mphi^2])/(32 \[Pi]^3 T)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


D\[CapitalGamma]\[Phi]to2chiEner[T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=D\[CapitalGamma]\[Phi]to2chiEnergy[T,mphi,mchi]


D\[CapitalGamma]\[Phi]to2chi=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(coupl Sqrt[1-(4 mchi^2)/mphi^2] mphi^2 BesselK[2,mphi/T] HeavisideTheta[1-(4 mchi^2)/mphi^2])/(16 \[Pi]^3 T)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


D\[CapitalGamma]\[Phi]to2\[Chi][T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=D\[CapitalGamma]\[Phi]to2chi[T,mphi,mchi]


Clear[cBEwithA];
cBEwithA[ms_,ls_,lhs_,\[Xi]inf_,mchi_?NumericQ,\[Kappa]_?NumericQ,ai_,af_(*,Ni_*)]:=
Module[{solStandard,aFstandard,TsStand,NsS,Yf,datTemps,datNes,TsstandardAsFuncOfTds,TsstandardAsFuncOfNds,xff,bgrcosmo},
(*Note the factor 1/(1+\[Alpha]) in the temperature evolution. One has to include this term if one wants to account for thermal corrections-> <(n/(2E))>((dm^2)/dt)*)
(*bgrcosmo=BgrCosmo[3/8,0,150.,0.3];
Tstandard=bgrcosmo[[1]];*)
solStandard=NDSolve[{Nn'[a]==(ls^3 averagedCif[ms/Ts[a],ms]/(a^7(*H[Tstandard[a]]*)H[Tstandard[a]]) Nn[a]^2 (Neq[Ts[a],ms,a]-Nn[a])+(lhs^2 a^2)/H[Tstandard[a]] (\[CapitalGamma]htophiphi[Tstandard[a],ms])-(\[Kappa]^2 a^2)/H[Tstandard[a]] Nn[a]/Neq[Ts[a],ms,a] \[CapitalGamma]\[Phi]to2\[Chi][Ts[a],ms,mchi]),
Ts'[a]==(lhs^2 /(a H[Tstandard[a]]) \[CapitalGamma]htophiphii2Energydens[Tstandard[a],ms]-\[Kappa]^2/(a H[Tstandard[a]]) Nn[a]/Neq[Ts[a],ms,a] \[CapitalGamma]\[Phi]to2chiEner[Ts[a],ms,mchi]-3/a ( Nn[a]Ts[a]/a^3)-Nn'[a]/Neq[Ts[a],ms,a] \[Rho]eq[Ts[a],ms])/(Nn[a]/Neq[Ts[a],ms,a] D\[Rho]eq[Ts[a],ms]-a^3 (\[Rho]eq[Ts[a],ms]/(Ts[a]Neq[Ts[a],ms,a]))^2 Nn[a]),Nn[ai]==Neq[\[Xi]inf Tstandard[ai],ms,ai](*Ni*),Ts[ai]==\[Xi]inf Tstandard[ai]
},{Nn,Ts},{a,ai,af//N},PrecisionGoal->9,AccuracyGoal->14,MaxSteps->10^7(*,Method->{"DiscontinuityProcessing"->False}*),MaxSteps->10^7];
aFstandard = solStandard[[1]][[1]][[2]]["Domain"][[1]][[2]];
TsStand[a_]:=Ts[a]/.solStandard[[1]];
NsS[a_]:=Nn[a]/.solStandard[[1]];
datTemps=Table[{ms/Tstandard[10^a],ms/TsStand[10^a]},{a,Log10[1],Log10[aFstandard],(Log10[aFstandard]-Log10[1])/1000}];
datNes=Table[{ms/Tstandard[10^a],NsS[10^a]/(10^(3a) s[Tstandard[10^a]])},{a,Log10[1],Log10[aFstandard],(Log10[aFstandard]-Log10[1])/1000}];
xff=datTemps[[-1]][[1]];
TsstandardAsFuncOfTds=Interpolation[datTemps,InterpolationOrder->1,Method->"Spline"];
TsstandardAsFuncOfNds=Interpolation[datNes,InterpolationOrder->1,Method->"Spline"];
Return[{TsstandardAsFuncOfTds(*=(x,xphi)*),TsstandardAsFuncOfNds(*(x,Yphi)*),TsStand(*=(a,Tphi)*),NsS(*=(a,Nphi)*),aFstandard,xff}]
];


Newton2[fi_,gF_,gJ_,{tol_,nmax_},flag_,x_]:=Block[{error=1,step=0,f=fi,v,rstep=0,t={},it,fit,errorit},

  While[error>tol && step++<nmax,
(*	v=Check[LinearSolve[gJ[f],gF[f]],"err"];If[v=="err",Abort[]];*)
	v=LinearSolve[gJ[f],gF[f]];
	f=(f-v);
	If[Min[Re[f]]<0,error=0.;Break[]];
	error=Max[Abs[gF[f]/f]];
	
	AppendTo[t,{error,f}];
	t=DeleteDuplicates[t,#1[[1]]==#2[[1]]&];
	If[Length[t]>=3,
	 it=Check[Interpolation[t,InterpolationOrder->1],"err"];
	 If[it!="err",
	 Quiet[fit=it[0.]];
	 If[Min[fit]>0,
	  errorit=Max[Abs[gF[fit]/fit]];
	  If[errorit<error,f=fit];
	 ];];
	];
   ];
 If[error>tol,f=0.fi; (*Print["Warning: f=0 returned. Newton inaccurate... ",{error,step}]*)(*AppendTo[n1,{error,step,x,fi,errorintt}]*)];
 If[Im[Total[f]]!=0 || Min[Re[f]]<0, f=0.fi; (*Print["Warning: f=0 returned. f negative or complex... "]*)(*AppendTo[n2,{error,step,x,fi,errorintt}]*)];
 f
];
(* clear internal names *)
(*Remove[fi,gF,gJ,tol,nmax,flag,x,error,step,f,v,rstep,t,it,fit,errorit];*)


(* simple ODE solver where inputs are RHS function F and Jacobian function J: AM2 midpoint *)
ODESolver2[F_,J_,finit_,{xmin_,xmax_},{hinit_,tol_,safety_}]:=Block[{xi,xip1,xf,fi,Fi,h,IM,gE,gT,gJ,fip1E,fip1T,error,h2,xip12,Fip12,flag},

xi=xmin;
xf=xmax;
fi=finit;
Fi=F[xi,fi];
h=hinit;
IM=IdentityMatrix[Length[fi]];

If[OPTION["Monitor"],r={}];
Monitor[
While[xi<xf,
xip1=xi+h;
If[xip1==xi,Print["Stepsize underflow! at x = "xi," Aborting..."];Abort[]];
If[xip1>xf,xip1=xf;h=xf-xi];

(* implicit Euler step *)
h2=h/2;
xip12=xi+h2;
gE[fip12_]:=fip12-fi-h2/2*(Fi+F[xip12,fip12]);
gJ[fip12_]:=IM-h2/2*J[xip12,fip12];
fip1E=Newton2[fi,gE,gJ,{tol/10,10},flag,xi];
(*If[Total[fip1E]==0,AppendTo[e1,{xi,error}]];*)

If[Total[fip1E]!=0,
 Fip12=F[xip12,fip1E];
 gE[fip1_]:=fip1-fip1E-h2/2*(Fip12+F[xip1,fip1]);
 gJ[fip1_]:=IM-h2/2*J[xip1,fip1];
 fip1E=Newton2[(*fi*)fip1E,gE,gJ,{tol/10,10},flag,xi];
(* If[Total[fip1E]==0,AppendTo[e2,{xi,error}]];*)

If[Total[fip1E]!=0,
  (* implicit Trapezoidal step *)
  gT[fip1_]:=fip1-fi-h/2*(Fi+F[xip1,fip1]);
  gJ[fip1_]:=IM-h/2*J[xip1,fip1];
  fip1T=Newton2[(*fi*)fip1E,gT,gJ,{tol/10,10},flag,xi];
(*  If[Total[fip1T]==0,AppendTo[e3,{xi,error}]];*)
 ];
];

(* relative error *)
If[Total[fip1T]!=0 && Total[fip1E]!=0,
 error=Max[Abs[(fip1T-fip1E)/fip1T]]/tol;
,
 error=100;
];
(* only if error<1 the step was successful *) 
 If[error==100,h=10^-1*h;Continue[]];  
 If[error>1,h=Max[safety/Sqrt[error],0.1]*h;(*AppendTo[e4,{xi,error}];*)Continue[]];  
 If[error==0.,h=10.*h,h=Min[safety/Sqrt[error],10.]*h];

(* step successful *)
xi = xip1;
fi=fip1E;
Fi=F[xi,fi];

(* if one of the particles effecticely decayed stop & restart as with only 1 state remaining *)
(*If[Min[Re[fi]]<SETTING["Y_min_value"] && xi>1,Break[]];*)
If[Re[fi[[1]]]<SETTING["Y_min_value"] && xi>1, FLAG["1 decayed"]=True;Break[]];
If[Length[fi]>1,
 If[Re[fi[[2]]]<SETTING["Y_min_value"] && xi>1, FLAG["2 decayed"]=True;Break[]];
];

Sow[Flatten[{xi,fi}]];
If[Head[r]==List,AppendTo[r,Flatten[{xi,fi}]]];
];
,
If[OPTION["Monitor"],
	{ListLogLogPlot[{r[[All,{1,2}]]//Re,r[[All,{1,3}]]//Re},ImageSize->Medium],
	 If[Length[r[[-1]]]==5,ListLogLogPlot[{r[[All,{1,4}]]//Re,r[[All,{1,5}]]//Re},ImageSize->Medium],Nothing]},
	{xi,fi}
]
];
];
(* clear internal names *)
(*Remove[F,J,finit,xmin,xmax,hinit,tol,safety,i,xip1,xf,fi,Fi,h,IM,gE,gT,gJ,fip1E,fip1T,error,h2,xip12,Fip12,flag];
*)


InterThermalAverages[];


z[Ns_,Ts_,a_,ms_]:=Ns/Neq[Ts,ms,a];

Clear[CE]
CE[Ns_,Ts_,a_,ms_,mchi_,ls_,\[Kappa]_,lhs_]:=lhs^2 \[CapitalGamma]htophiphii2Energydens[Tstandard[a],ms]-\[Kappa]^2 z[Ns,Ts,a,ms] \[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi]


DzN[Ts_,a_,ms_]:=1/Neq[Ts,ms,a];
DzTs[Ns_,Ts_,a_,ms_]:=-Ns \[Rho]eq[Ts,ms]/(Neq[Ts,ms,a]^2 Ts^2) a^3;


DaveragedCif[xs_,ms_]:=Module[{\[Epsilon]= xs/1000},(averagedCif[xs+\[Epsilon],ms]-averagedCif[xs-\[Epsilon],ms])/(2\[Epsilon])]


DaveragedCT[Ts_,ms_]:=(-ms/Ts^2)*DaveragedCif[ms/Ts,ms]


Clear[F]
F[Ns_,Ts_,a_,ms_,mchi_,ls_,\[Kappa]_,lhs_]:=Block[{Tsprime,Nprime},
Nprime =ls^3 averagedCif[ms/Ts,ms]/(a^7 H[Tstandard[a]]) Ns^2 (Neq[Ts,ms,a]-Ns)+a^2/H[Tstandard[a]] (lhs^2 (\[CapitalGamma]htophiphi[Tstandard[a],ms])-\[Kappa]^2 Ns/Neq[Ts,ms,a] \[CapitalGamma]\[Phi]to2\[Chi][Ts,ms,mchi]);
Tsprime =(CE[Ns,Ts,a,ms,mchi,ls,\[Kappa],lhs]/(a H[Tstandard[a]])-3 Ns Ts/a^4-Nprime/Neq[Ts,ms,a] \[Rho]eq[Ts,ms])/(z[Ns,Ts,a,ms]D\[Rho]eq[Ts,ms]-a^3 (\[Rho]eq[Ts,ms]/(Ts Neq[Ts,ms,a]))^2 Ns);
Return[{Nprime,Tsprime}]
]


Clear[J]
J[Ns_,Ts_,a_,ms_,mchi_,ls_,\[Kappa]_,lhs_]:=Block[{NprimeN,NprimeT,TprimeN,TprimeT, Nprime},
Nprime=ls^3 averagedCif[ms/Ts,ms]/(a^7 H[Tstandard[a]]) Ns^2 (Neq[Ts,ms,a]-Ns)+a^2/H[Tstandard[a]] (lhs^2 (\[CapitalGamma]htophiphi[Tstandard[a],ms])-\[Kappa]^2 Ns/Neq[Ts,ms,a] \[CapitalGamma]\[Phi]to2\[Chi][Ts,ms,mchi]);

(*********************************************************************************)
(*********************************************************************************)
NprimeN=-((ls^3 Ns^2 averagedCif[ms/Ts,ms])/(a^7 H[Tstandard[a]]))+(2 ls^3 Ns averagedCif[ms/Ts,ms] (-Ns+Neq[Ts,ms,a]))/(a^7 H[Tstandard[a]])-(a^2 \[Kappa]^2 \[CapitalGamma]\[Phi]to2\[Chi][Ts,ms,mchi])/(H[Tstandard[a]] Neq[Ts,ms,a]);
NprimeT = (ls^3 Ns^2 DaveragedCT[Ts,ms] (-Ns+Neq[Ts,ms,a]))/(a^7 H[Tstandard[a]])+(ls^3 Ns^2 averagedCif[ms/Ts,ms] \[Rho]eq[Ts,ms])/(a^4 Ts^2 H[Tstandard[a]])+(a^2 (-((Ns \[Kappa]^2 D\[CapitalGamma]\[Phi]to2\[Chi][Ts,ms,mchi])/Neq[Ts,ms,a])+(a^3 Ns \[Kappa]^2 \[CapitalGamma]\[Phi]to2\[Chi][Ts,ms,mchi] \[Rho]eq[Ts,ms])/(Ts^2 Neq[Ts,ms,a]^2)))/H[Tstandard[a]];
TprimeN=-(((-((3 Ns Ts)/a^4)+1/(a H[Tstandard[a]]) (lhs^2 \[CapitalGamma]htophiphii2Energydens[Tstandard[a],ms]-\[Kappa]^2 z[Ns,Ts,a,ms] \[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi])-(Nprime \[Rho]eq[Ts,ms])/Neq[Ts,ms,a]) (D\[Rho]eq[Ts,ms]/Neq[Ts,ms,a]-(a^3 \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2)))/(D\[Rho]eq[Ts,ms] z[Ns,Ts,a,ms]-(a^3 Ns \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2))^2)+(-((3 Ts)/a^4)-(\[Kappa]^2 \[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi])/(a H[Tstandard[a]] Neq[Ts,ms,a])-(NprimeN \[Rho]eq[Ts,ms])/Neq[Ts,ms,a])/(D\[Rho]eq[Ts,ms] z[Ns,Ts,a,ms]-(a^3 Ns \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2));
TprimeT=-(((-((3 Ns Ts)/a^4)+1/(a H[Tstandard[a]]) (lhs^2 \[CapitalGamma]htophiphii2Energydens[Tstandard[a],ms]-\[Kappa]^2 z[Ns,Ts,a,ms] \[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi])-(Nprime \[Rho]eq[Ts,ms])/Neq[Ts,ms,a]) (DD\[Rho]eq[Ts,ms] z[Ns,Ts,a,ms]-(3 a^3 Ns D\[Rho]eq[Ts,ms] \[Rho]eq[Ts,ms])/(Ts^2 Neq[Ts,ms,a]^2)+(2 a^3 Ns \[Rho]eq[Ts,ms]^2)/(Ts^3 Neq[Ts,ms,a]^2)+(2 a^6 Ns \[Rho]eq[Ts,ms]^3)/(Ts^4 Neq[Ts,ms,a]^3)))/(D\[Rho]eq[Ts,ms] z[Ns,Ts,a,ms]-(a^3 Ns \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2))^2)+(-((3 Ns)/a^4)-(Nprime D\[Rho]eq[Ts,ms])/Neq[Ts,ms,a]-(NprimeT \[Rho]eq[Ts,ms])/Neq[Ts,ms,a]+(a^3  \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2) Nprime+(-\[Kappa]^2 D\[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi] z[Ns,Ts,a,ms]+(a^3 Ns \[Kappa]^2 \[CapitalGamma]\[Phi]to2chiEner[Ts,ms,mchi] \[Rho]eq[Ts,ms])/(Ts^2 Neq[Ts,ms,a]^2))/(a H[Tstandard[a]]))/(D\[Rho]eq[Ts,ms] z[Ns,Ts,a,ms]-(a^3 Ns \[Rho]eq[Ts,ms]^2)/(Ts^2 Neq[Ts,ms,a]^2));
Return[Re[{{NprimeN,NprimeT},{TprimeN,TprimeT}}]]

]


Clear[JacNum]
JacNum[Ns_?NumericQ,Ts_?NumericQ,a_?NumericQ,ms_,mchi_,ls_,\[Kappa]_,lhs_]:=Module[{hN,hT},hN=10^-7 Max[Abs[Ns],1.];
hT=10^-7 Max[Abs[Ts],1.];
Re[{(F[Ns+hN,Ts,a,ms,mchi,ls,\[Kappa],lhs]-F[Ns-hN,Ts,a,ms,mchi,ls,\[Kappa],lhs])/(2 hN),(F[Ns,Ts+hT,a,ms,mchi,ls,\[Kappa],lhs]-F[Ns,Ts-hT,a,ms,mchi,ls,\[Kappa],lhs])/(2 hT)}\[Transpose]]];


afFromXf[ms_?NumericQ,xf_?NumericQ,ai_?NumericQ]:=Module[{u},u/. FindRoot[ms/Tstandard[Exp[u]]-xf==0,{u,Log[ai]+5},PrecisionGoal->12,AccuracyGoal->12]//Exp]


Clear[cBEsolver];
cBEsolver[mphi_?NumericQ,lphi_?NumericQ,lhphi_?NumericQ,\[Xi]inf_?NumericQ,mchi_?NumericQ,\[Kappa]_?NumericQ,xf_?NumericQ]:=
Module[{init,FF,JJ,FFu,JJu,u0,uf,r,rA,datTs,datNs,datY,TsOfx,NsOfx,YOfx,xff,af,ai=1,xi, rat},
xi=mphi/Tstandard[ai];
af=afFromXf[mphi,xf,ai];init={Neq[\[Xi]inf Tstandard[ai],mphi,ai],\[Xi]inf Tstandard[ai]};
FF[a_,Y_]:=F[Y[[1]],Y[[2]],a,mphi,mchi,lphi,\[Kappa],lhphi];
JJ[a_,Y_]:=J[Y[[1]],Y[[2]],a,mphi,mchi,lphi,\[Kappa],lhphi];
FFu[u_,Y_]:=Module[{a=Exp[u]},a*FF[a,Y]];
JJu[u_,Y_]:=Module[{a=Exp[u]},a*JJ[a,Y]];
u0=Log[ai];
uf=Log[af];
r=Reap[ODESolver2[FFu,JJu,init,{u0,uf},{0.05,0.001,0.9}]][[2,1]];
rA=r/. {u_?NumericQ,Ns_?NumericQ,Ts_?NumericQ}:>Module[{a=Exp[u],Tsm=Tstandard[Exp[u]]},{mphi/Tsm,a,Ns,Ts}];
datTs=SortBy[rA[[All,{1,4}]],First];(*{x,Ts}*)
datNs=SortBy[rA[[All,{1,3}]],First];(*{x,Ns}*)
datY=SortBy[rA/. {x_,a_,Ns_,Ts_}:>{x,(Ns/a^3)/s[Tstandard[a]]},First];
xff=datTs[[-1,1]];
TsOfx=Interpolation[datTs,InterpolationOrder->1];
NsOfx=Interpolation[datNs,InterpolationOrder->1];
YOfx=Interpolation[datY,InterpolationOrder->1];
rat= datY[[-10]][[2]]/Yobs[mphi];
Return[{TsOfx,NsOfx,YOfx,{xi,xf},rat}]]



Clear[cBEsolverwithA];
cBEsolverwithA[mphi_?NumericQ,lphi_?NumericQ,lhphi_?NumericQ,\[Xi]inf_?NumericQ,mchi_?NumericQ,\[Kappa]_?NumericQ,xf_?NumericQ]:=
Module[{init,FF,JJ,FFu,JJu,u0,uf,r,rA,datTs,datNs,datY,TsOfA,NsOfA,YOfA,xff,af,ai=1,xi},
xi=mphi/Tstandard[ai];
af=afFromXf[mphi,xf,ai];init={Neq[\[Xi]inf Tstandard[ai],mphi,ai],\[Xi]inf Tstandard[ai]};
FF[a_,Y_]:=F[Y[[1]],Y[[2]],a,mphi,mchi,lphi,\[Kappa],lhphi];
JJ[a_,Y_]:=J[Y[[1]],Y[[2]],a,mphi,mchi,lphi,\[Kappa],lhphi];
FFu[u_,Y_]:=Module[{a=Exp[u]},a*FF[a,Y]];
JJu[u_,Y_]:=Module[{a=Exp[u]},a*JJ[a,Y]];
u0=Log[ai];
uf=Log[af];
r=Reap[ODESolver2[FFu,JJu,init,{u0,uf},{0.05,0.001,0.9}]][[2,1]];
rA=r/. {u_?NumericQ,Ns_?NumericQ,Ts_?NumericQ}:>Module[{a=Exp[u],Tsm=Tstandard[Exp[u]]},{mphi/Tsm,a,Ns,Ts}];
datTs=SortBy[rA[[All,{2,4}]],First];(*{a,Ts}*)
datNs=SortBy[rA[[All,{2,3}]],First];(*{a,Ns}*)
datY=SortBy[rA/. {x_,a_,Ns_,Ts_}:>{a,(Ns/a^3)/s[Tstandard[a]]},First];
xff=datTs[[-1,1]];
TsOfA=Interpolation[datTs,InterpolationOrder->1];
NsOfA=Interpolation[datNs,InterpolationOrder->1];
YOfA=Interpolation[datY,InterpolationOrder->1];
Return[{TsOfA,NsOfA,YOfA,{xi,xf}}];]


Clear[nBEwithA];
nBEwithA[ms_,\[Xi]inf_,mchi_?NumericQ,\[Kappa]_?NumericQ,Tphi_,zphi_,ai_,af_(*,Ni_*)]:=
Module[{solStandard,aFstandard,NsS,Yf,datTemps,datYes,TsstandardAsFuncOfTds,TsstandardAsFuncOfNds,xff,rat},
solStandard=NDSolve[{Nn'[a]==((\[Kappa]^2 a^2)/H[Tstandard[a]] zphi[a] \[CapitalGamma]\[Phi]to2\[Chi][Tphi[a],ms,mchi]),Nn[ai]==Neq[\[Xi]inf Tstandard[ai],ms,ai]
},Nn,{a,ai,af//N},PrecisionGoal->9,AccuracyGoal->14,MaxSteps->10^7,MaxSteps->10^7];
aFstandard = solStandard[[1]][[1]][[2]]["Domain"][[1]][[2]];
NsS[a_]:=Nn[a]/.solStandard[[1]];
datYes=Table[{ms/Tstandard[10^a],NsS[10^a]/(10^(3a) s[Tstandard[10^a]])},{a,Log10[1],Log10[aFstandard],(Log10[aFstandard]-Log10[1])/1000}];
xff=datYes[[-1]][[1]];
TsstandardAsFuncOfNds=Interpolation[datYes,InterpolationOrder->1,Method->"Spline"];
rat=2 datYes[[-30]][[2]]/Yobs[mchi];
Return[{TsstandardAsFuncOfNds(*(x,Yphi)*),NsS(*=(a,Nphi)*),aFstandard,xff,datYes[[-30]][[2]],rat}]
];


Clear[RatedecayC0]
RatedecayC0[Ts_,Y_,k_,mphi_,mchi_]:=Module[{ratec0},
ratec0[x_]:=k^2/(neq[Ts[x],mphi]) \[CapitalGamma]\[Phi]to2\[Chi][Ts[x],mphi,mchi];
Return[ratec0]
];


Clear[RatedecayC2]
RatedecayC2[Ts_,Y_,k_,mphi_,mchi_]:=Module[{ratec2},
ratec2[x_]:=k^2/(neq[Ts[x],mphi]) Ts[x]/(3mphi) \[CapitalGamma]\[Phi]to2\[Chi]C2[Ts[x],mphi,mchi];
Return[ratec2]
]


(* ::Text:: *)
(*Cannibal*)


(*Clear[RateCannibal3to2]*)
RateCannibal3to2[Ts_,Y_,ls_,mphi_]:=Module[{rate3to2},
rate3to2[x_]:=ls^3 averagedCif[mphi/Ts[x],mphi]s[mphi/x]^2 Y[x]^2;
Return[rate3to2]
];


(*Clear[RateCannibal2to3]*)
RateCannibal2to3[Ts_,Y_,ls_,mphi_]:=Module[{rate2to3},
rate2to3[x_]:=ls^3 averagedCif[mphi/Ts[x],mphi]s[mphi/x]Y[x]neq[Ts[x],mphi] ;
Return[rate2to3]
]


Delta[mphi_, mchi_] := Sqrt[mphi^2 - 4 mchi^2];
Msq[mphi_, mchi_, y_] := 2 y^2 (mphi^2 - 4 mchi^2);

pphys[afinal_, p_] := (afinal) p;
Echi[mchi_,q_,a_] := Sqrt[mchi^2 + (q/a)^2];

Ephimin[mphi_, mchi_,q_,a_] := 
  (mphi^2 Echi[mchi,q,a] - 
    mphi (q/a) Delta[mphi, mchi])/(2 mchi^2);

Ephimax[mphi_, mchi_,q_,a_] := 
  (mphi^2 Echi[mchi,q,a] + 
    mphi (q/a) Delta[mphi, mchi])/(2 mchi^2);

Kernel[mphi_, mchi_, y_,Tphi_,z_,q_,a_] := z[a](Msq[mphi, mchi, y]/(8 Pi))*
  (Tphi[a]/(Echi[mchi,q,a] (q/a)))* (Exp[-Ephimin[mphi, mchi,q,a]/Tphi[a]] - Exp[-Ephimax[mphi, mchi,q,a]/Tphi[a]]);

fc[mphi_, mchi_, y_,Tphi_,z_,q_,a_] := NIntegrate[
  Kernel[mphi, mchi, y,Tphi,z,q,ap]/(ap H[Tstandard[ap]]),
  {ap, 1, a}(*, WorkingPrecision->MachinePrecision,PrecisionGoal->15,AccuracyGoal->15*)
];


fchi[mphi_, mchi_, y_,Tphi_,z_,p_,a_]:= fc[mphi, mchi, y,Tphi,z,a p,a] 
