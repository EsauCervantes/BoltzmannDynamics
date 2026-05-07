(* ::Package:: *)

(* ::Section:: *)
(*Definition of collision integral routines*)


(* ::Text:: *)
(*Here I defined the matrix element in the Lab frame, p1,p2->p3,p4,p5 with p3+p4 at rest. It turns out that the integrals are much faster this way than in the CM frame (p1+p2 at rest).*)


SetDirectory[NotebookDirectory[]];


M=Compile[{{st,_Real,0},{E5t,_Real,0},{x,_Real,0},{xs,_Real,0},{\[Phi],_Real,0}},
Block[{p3t,p4t,p5t,p1t,p2t,t13p,t14p,t15p,s35p,s45p,s34p,t25p,t23p,t24p},(*Here everything is already boosted from the lab frame (CM of two final states) to the total CM frame*)
p3t={-(((-3+E5t Sqrt[st]) (9-6 E5t Sqrt[st]+st+Sqrt[3] Sqrt[((-1+E5t^2) st)/(-3+E5t Sqrt[st])^2] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] xs))/(2 Sqrt[st] (9-6 E5t Sqrt[st]+st))),(Sqrt[-1+E5t^2] (-3+E5t Sqrt[st]) Sqrt[1-x^2])/(2 Abs[-3+E5t Sqrt[st]])+(Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] (E5t st Sqrt[1-x^2] xs-3 Sqrt[st-st x^2] xs-Sqrt[st (9-6 E5t Sqrt[st]+st)] x Sqrt[1-xs^2] Cos[\[Phi]]))/(2 st (9-6 E5t Sqrt[st]+st)),(Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] Sqrt[1-xs^2] Sin[\[Phi]])/(2 Sqrt[st (9-6 E5t Sqrt[st]+st)]),(Sqrt[-1+E5t^2] (-3+E5t Sqrt[st]) x)/(2 Abs[-3+E5t Sqrt[st]])+(Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] (-3 Sqrt[st] x xs+E5t st x xs+Sqrt[st (9-6 E5t Sqrt[st]+st)] Sqrt[1-x^2] Sqrt[1-xs^2] Cos[\[Phi]]))/(2 st (9-6 E5t Sqrt[st]+st))};
p4t={-(1/(2 (9-6 E5t Sqrt[st]+st)))(E5t-3/Sqrt[st]) (9-6 E5t Sqrt[st]+st-Sqrt[3] Sqrt[((-1+E5t^2) st)/(-3+E5t Sqrt[st])^2] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] xs),(Sqrt[-1+E5t^2] (-3+E5t Sqrt[st]) Sqrt[1-x^2])/(2 Abs[-3+E5t Sqrt[st]])+(Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] (-E5t st Sqrt[1-x^2] xs+3 Sqrt[st-st x^2] xs+Sqrt[st (9-6 E5t Sqrt[st]+st)] x Sqrt[1-xs^2] Cos[\[Phi]]))/(2 st (9-6 E5t Sqrt[st]+st)),-((Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] Sqrt[1-xs^2] Sin[\[Phi]])/(2 Sqrt[st (9-6 E5t Sqrt[st]+st)])),(Sqrt[-1+E5t^2] (-3+E5t Sqrt[st]) x)/(2 Abs[-3+E5t Sqrt[st]])-(Sqrt[3] Sqrt[27+4 E5t (-9+st) Sqrt[st]-6 st+12 E5t^2 st-st^2] (-3 Sqrt[st] x xs+E5t st x xs+Sqrt[st (9-6 E5t Sqrt[st]+st)] Sqrt[1-x^2] Sqrt[1-xs^2] Cos[\[Phi]]))/(2 st (9-6 E5t Sqrt[st]+st))};
p5t={E5t,Sqrt[-1+E5t^2] Sqrt[1-x^2],0,Sqrt[-1+E5t^2] x};
p1t={3/(2 Sqrt[st]),0,0,1/2 Sqrt[-4+9/st]};
p2t={3/(2 Sqrt[st]),0,0,-(1/2) Sqrt[-4+9/st]};
t13p=2-2(p1t[[1]]p3t[[1]]-p1t[[2]]p3t[[2]]-p1t[[3]]p3t[[3]]-p1t[[4]]p3t[[4]]);
t14p=2-2(p1t[[1]]p4t[[1]]-p1t[[2]]p4t[[2]]-p1t[[3]]p4t[[3]]-p1t[[4]]p4t[[4]]);
t15p=2-2(p1t[[1]]p5t[[1]]-p1t[[2]]p5t[[2]]-p1t[[3]]p5t[[3]]-p1t[[4]]p5t[[4]]);
s35p=2+2(p3t[[1]]p5t[[1]]-p3t[[2]]p5t[[2]]-p3t[[3]]p5t[[3]]-p3t[[4]]p5t[[4]]);
s45p=2+2(p4t[[1]]p5t[[1]]-p4t[[2]]p5t[[2]]-p4t[[3]]p5t[[3]]-p4t[[4]]p5t[[4]]);
s34p=(9-6 E5t Sqrt[st]+st)/st;
t25p=2-2(p2t[[1]]p5t[[1]]-p2t[[2]]p5t[[2]]-p2t[[3]]p5t[[3]]-p2t[[4]]p5t[[4]]);
t23p=2-2(p2t[[1]]p3t[[1]]-p2t[[2]]p3t[[2]]-p2t[[3]]p3t[[3]]-p2t[[4]]p3t[[4]]);
t24p=2-2(p2t[[1]]p4t[[1]]-p2t[[2]]p4t[[2]]-p2t[[3]]p4t[[3]]-p2t[[4]]p4t[[4]]);
(*Sqrt[3] /m*) (1/(1-s34p)+1/(1-s35p)+1/(1-s45p)+st/(-9+st)+(3 (3+s35p (-2+s45p)-2 s45p+s34p (-2+s35p+s45p)) st)/((-1+s34p) (-1+s35p) (-1+s45p) (-9+st))+1/(1-t13p)+1/(1-t14p)-3 (1/((-1+s45p) (-1+t13p))+1/((-1+s35p) (-1+t14p))+1/((-1+s34p) (-1+t15p)))+1/(1-t15p)+1/(1-t23p)+1/(1-t24p)-3 (1/((-1+s45p) (-1+t23p))+1/((-1+s35p) (-1+t24p))+1/((-1+s34p) (-1+t25p)))-3 (1/((-1+t14p) (-1+t23p))+1/((-1+t15p) (-1+t23p))+1/((-1+t13p) (-1+t24p))+1/((-1+t15p) (-1+t24p))+(-2+t13p+t14p)/((-1+t13p) (-1+t14p) (-1+t25p)))+1/(1-t25p))
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


Integrand=With[{MM=M},
Compile[{{st,_Real,0},{E5t,_Real,0},{x,_Real,0},{xs,_Real,0},{\[Phi],_Real,0}},
(*3/m^2*)MM[st,E5t,x,xs,\[Phi]]^2(*Sqrt[3]m*)( Sqrt[-1+E5t^2] Sqrt[3-2 E5t Sqrt[st]-st])/Sqrt[9-6 E5t Sqrt[st]+st](*Overall a factor of Sqrt[3]3/m*),
RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
]
];
\[Sigma]2to3[m_?NumericQ,s_?NumericQ,points_,ac_]:=1/2*1/(2Sqrt[-4 m^2 s+s^2])*1/(8(2Pi)^4) 3/m^2 Sqrt[3]m m NIntegrate[Integrand[9 m^2/s,E5t,x,xs,\[Phi]],{E5t,1,1/(2m) (-3 m^2+s)/Sqrt[s]},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},WorkingPrecision->MachinePrecision,AccuracyGoal->ac,PrecisionGoal->15(*,
Method\[Rule]{"QuasiMonteCarlo","MaxPoints"->points}*)](*1/3!Symmetry factors*)
nphieq=Compile[{{x,_Real,0},{m,_Real,0}},
1/(x*2Pi^2)*m^3*BesselK[2,x],
RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];
facC0=1/(2(2Pi)^4);

average2to3[xt_,m_,ac_,points_]:=1/(8(2Pi)^4)*1/4 facC0 *m/xt 9m^2 m *(Sqrt[3]3)/m*m *1/nphieq[xt,m]^3 NIntegrate[1/(*Sqrt[-4 m^2 +s]*)1 1/st^2*BesselK[1,(3xt)/Sqrt[st]]*(*Sqrt[s]*)Sqrt[9 1/st-4]Integrand[st,E5t,x,xs,\[Phi]],{st,0,1/10^5,1/10^4,1/10^3,1/10^2,1/10,0.5,0.9,1},{E5t,1,-((-3+st)/(2 Sqrt[st]))},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},WorkingPrecision->MachinePrecision,AccuracyGoal->ac,PrecisionGoal->15,
Method->{"QuasiMonteCarlo","MaxPoints"->points}]


\[Sigma]2to3[0.1,9*0.1^2*1.5,10,10]


9*0.1^2


Ecm^2->s


Sqrt[9*0.1^2*1.2]/2


E1=Sqrt[9*0.1^2*1.5]/2;


4Pi 1/(64Pi^2 s)/.s->(2E1)^2


(0.0001920216100530995`/0.14736568804805122`)^-1


(*f\[Sigma][]:=Module[{lphi=1,mphi=0.1,g=Sqrt[3 *0.1^2],m,b,r},

];
f\[Sigma][];*)
lphi=1;mphi=0.1;g=Sqrt[3 *0.1^2];
r=Import["files/tabulated_sigma,lphi=1,mphi=0.1.dat"];
\[Sigma]i=Interpolation[r/.{x_,y_}->{x,Log10[y]},Method->"Spline",InterpolationOrder->1];
m=(\[Sigma]i[(35*mphi)^2]-Log10[\[Sigma]2to3[mphi,(70mphi)^2,6*10^6,19]])/(Log10[(35*mphi)^2]-Log10[(70*mphi)^2]);
b=\[Sigma]i[(35*mphi)^2]-m*Log10[(35*mphi)^2];
\[Sigma]f[s_]:=Piecewise[{{10^\[Sigma]i[s],s<=(35*0.1)^2},{10^(m Log10[s]+b),s>(35*0.1)^2}}]
\[Sigma]ff[st_?NumericQ]:=(*(0.1/mphi)^2*)\[Sigma]f[(9(0.1)^2)/st];


IntegrandAveraged2to3i=With[{sf=\[Sigma]ff},
Compile[{{st,_Real,0},{x,_Real,0}},
(*3 m^3*) (9-4 st) BesselK[1,(3 x)/Sqrt[st]] sf[st]/st^(3/2)
,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
]
];
average2to3i[x_?NumericQ]:=(*facC01/nphieq[x,m]^3m/x9m^2(0.1/m)^23m^3 *)NIntegrate[1/st^2 IntegrandAveraged2to3i[st,x],{st,0,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.8,0.85,0.97,0.9,0.92,0.95,0.99,0.999,1},WorkingPrecision->MachinePrecision,AccuracyGoal->30,PrecisionGoal->15];


fEp2t=Compile[{{x,_Real,0},{st,_Real,0}},
(*3m^3*)(  (Sqrt[9-4 st] BesselK[2,(3 x)/Sqrt[st]])/(st x)+2/Sqrt[st] NIntegrate[(E^(-((3 x)/(Ept Sqrt[st]))) Log[(3-Sqrt[1-Ept^2] Sqrt[9-4 st])/(3+Sqrt[1-Ept^2] Sqrt[9-4 st])])/Ept^2,{Ept,0,10^-5,10^-4,10^-2,10^-1,0.2,0.5,0.8,0.85,0.9,0.99,1},WorkingPrecision->MachinePrecision,AccuracyGoal->15,PrecisionGoal->15,Method->{Automatic,"SymbolicProcessing"->0}(*,AccuracyGoal->10,PrecisionGoal\[Rule]15*)]),
RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];
IntegrandAveraged2to3wP2Ei=With[{fEp22=fEp2t,sf=\[Sigma]ff},
Compile[{{x,_Real,0},{st,_Real,0}},
(*3 m^3*)fEp22[x,st]*(*6m^2*)Sqrt[9-4 st]/st*1/st^2 sf[st],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
]
];
fac1=1/(8(2Pi)^4);
average2to3wP2Ei[x_?NumericQ]:=(*fac11/nphieq[x,m]^39m^2(0.1/m)^2(6m^2)(3m^3)*)NIntegrate[IntegrandAveraged2to3wP2Ei[x,st],{st,0,10^-5,10^-4,10^-2,10^-1,0.2,0.5,0.8,0.85,0.9,0.99,0.999,1},AccuracyGoal->25,PrecisionGoal->25,WorkingPrecision->MachinePrecision]


average2to3wP2Ei[1.]//Timing


LogLogPlot[Qii[x,0.1],{x,0,20}]


(* ::Section:: *)
(*definition of Int[M^2 p^2/E] (cross section with factor p^2/E inside the final phase space integral).*)


dataQ=Import["files/eta_integral2.dat"];
Qi= Interpolation[dataQ,Method->"Spline",InterpolationOrder->2];
Qii=Compile[{{a,_Real,0},{b,_Real,0}},Piecewise[{{10^((Log10[Qi[0.003,b]]-Log10[Qi[0.0035,b]])/(Log10[0.003]-Log10[0.0035]) Log10[a]+Log10[Qi[0.003,b]]-Log10[0.003] (Log10[Qi[0.003,b]]-Log10[Qi[0.0035,b]])/(Log10[0.003]-Log10[0.0035])),a<0.003},{Qi[a,b],0.003<=a<=22},{Re[10^((Log10[Qi[22,b]]-Log10[Qi[23,b]])/(22-23) a+Log10[Qi[22,b]]-22*(Log10[Qi[22,b]]-Log10[Qi[23,b]])/(22-23))],a>22}}],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];
IntegrandP2Averaged=With[{QQi=Qii,MM=M},
Compile[{{xT,_Real,0},{st,_Real,0},{E5t,_Real,0},{xA,_Real,0},{xs,_Real,0},{\[Phi],_Real,0}},

(*m^4/6*)(E^(-((3 xT)/Sqrt[st]))  Sqrt[-(4/3)+3/st] Sqrt[-((-3+2 E5t Sqrt[st]+st)/st)] (E5t (2 (-1+E5t^2) st xA (Sqrt[st]+3 xT)+9 E^((3 xT)/Sqrt[st]) E5t Sqrt[(-1+E5t^2) st] xT^2 BesselK[2,(3 xT)/Sqrt[st]])-27 E^((3 xT)/Sqrt[st]) Sqrt[-1+E5t^2] xT^3 QQi[(3 xT)/Sqrt[st],(Sqrt[-1+E5t^2] xA)/E5t]))/(E5t Sqrt[9-6 E5t Sqrt[st]+st] xT^3)(*3 /m^2*)MM[st,E5t,xA,xs,\[Phi]]^2 1/st^2(*Altogether we have a factor of 3/m^2 from MM^2 and m^4/6 from the other terms, hence the final result should be multiplied by 1/2m^2*)
,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True]
];
fac=1/8*(4Pi)^2/(2Pi)^6*1/(2^4 (2Pi)^4);
factor=(1/(2Pi)^4 1/4)*1/2*1/(8(2Pi)^4);
P2Eaveraged2[xT_?NumericQ,points_?NumericQ,ac_?NumericQ]:=(*1/nphieq[xT,m]^39m^2mm^2/2*)factor (*1/nphieq[xT,1]^3*)NIntegrate[(*IntegrandP2Averaged[xT,st,E5t,xA,xs,\[Phi]]*)(E^(-((3 xT)/Sqrt[st]))  Sqrt[-(4/3)+3/st] Sqrt[-((-3+2 E5t Sqrt[st]+st)/st)] (E5t (2 (-1+E5t^2) st xA (Sqrt[st]+3 xT)+9 E^((3 xT)/Sqrt[st]) E5t Sqrt[(-1+E5t^2) st] xT^2 BesselK[2,(3 xT)/Sqrt[st]])-27 E^((3 xT)/Sqrt[st]) Sqrt[-1+E5t^2] xT^3 Qii[(3 xT)/Sqrt[st],(Sqrt[-1+E5t^2] xA)/E5t]))/(E5t Sqrt[9-6 E5t Sqrt[st]+st] xT^3)(*3 /m^2*)M[st,E5t,xA,xs,\[Phi]]^2 1/st^2,{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,1},{E5t,1,-( (-3+st)/(2 Sqrt[st]))},{xA,-1,0,1},{xs,-1,0,1},{\[Phi],0,2Pi},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->13,WorkingPrecision->MachinePrecision];


st/Sqrt[9-4 st]*1/(1536 (m^2) (\[Pi]^4) )


factor=(1/(2Pi)^4 1/4)*1/2*1/(8(2Pi)^4);


factor


1/(16384 \[Pi]^8) (9 m^5)/2


\[Sigma]T[st_?NumericQ,cosh_?NumericQ,points_,ac_]:=(*the first 1/2 comes from integrating two momenta in their CM. The second term comes from (1/(4F)). The third term comes from 1/(2^3(2\[Pi])^9)(2\[Pi])(2\[Pi])^4*)st/ Sqrt[9-4 st](*1/(1536 m^2 \[Pi]^4)*Sqrt[3]m1/2lphi^3/m^2m m  <- one m from dE5 = m dE5t and one from p3^2/E3*) NIntegrate[Integrand[st,E5t,x,xs,\[Phi]](cosh E5t+Sqrt[-1+cosh^2] Sqrt[-1+E5t^2] x-1/(cosh E5t+Sqrt[-1+cosh^2] Sqrt[-1+E5t^2] x)),{E5t,1,(3-st)/(2 Sqrt[st])},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},WorkingPrecision->MachinePrecision,AccuracyGoal->ac,PrecisionGoal->20,
Method->{"QuasiMonteCarlo","MaxPoints"->points}]


f\[Sigma]T[]:=
Module[{r},
r={};
(*Table[{10^st//N,10^cosh//N,\[Sigma]T[k,10^st,10^cosh,10^6,15]},{st,Log10[1/150],Log10[0.99999],(Log10[0.99999]-Log10[1/150])/40},{cosh,Log10[1.],Log10[10],(Log10[10.]-Log10[1.])/40}]*)
Do[AppendTo[r,{10^-stt//N,10^cosh//N,\[Sigma]T[10^-stt,10^cosh,10^6,15]}];Export["files/tabulated_sigma,second_m.dat",r],{stt,Log10[1.000001],Log10[3],(Log10[3]-Log10[1.000001])/55},{cosh,Log10[1.],Log10[6],(Log10[6]-Log10[1])/35}]
(*r=Import["files/tabulated_sigmaT,k="<>ToString[k]<>",z3.dat"];
\[Sigma]i=Interpolation[r/.{x_,y_}->{x,Log10[y]},Method->"Spline",InterpolationOrder\[Rule]1];(*This trick is to avoid the interpolation to be negative for some points*)
m=(\[Sigma]i[1/150]-Log10[\[Sigma]2to3z3[1,1/1000,6*10^6,19]])/(Log10[1/150]-Log10[1/1000]);
b=\[Sigma]i[1/150]-m*Log10[1/150];
\[Sigma]f[st_]:=Piecewise[{{10^\[Sigma]i[st],st>1/150},{10^(m Log10[st]+b),st<1/150}}]*)
];


f\[Sigma]T[]


\[Sigma]tabulated=Import["files/tabulated_sigma,second_m.dat"];


Pij=Sqrt[s-4m^2]/2;


Assuming[m>0&&st>0&&E5t>0,st/Sqrt[9-4 st](*the term on the right is 1/(4F) in the cross section factor*) Pij *Sqrt[s]4 Pij *Sqrt[s]/.s->9 m^2/st//FullSimplify]


Assuming[m>0&&st>0,Exp[-Sqrt[s] cosh/Tphi]/.s->9 m^2/st/.Tphi->m/xT//FullSimplify]


\[Sigma]Ti = Interpolation[\[Sigma]tabulated/.{x_,y_,z_}->{x,y,Log10[z]},InterpolationOrder->1,Method->"Spline"];


1/(1536 m^2 \[Pi]^4) Sqrt[3]m 3 /m^2 m


1/(1536 m^2 \[Pi]^4)*Sqrt[3] 3/m m m


1/(1536 m^2 \[Pi]^4)*Sqrt[3] 3/m m m 1/(2(2Pi)^4) 9m^4 9m^2


averagesecondmom1[xT_?NumericQ,ac_,points_]:=NIntegrate[Sqrt[9-4 st]/st^3 Integrand[st,E5t,x,xs,\[Phi]]((Sqrt[st] BesselK[2,(3 xT)/Sqrt[st]])/(3 xT) E5t+2/27 (E^(-((3 xT)/Sqrt[st])) st (Sqrt[st]+3 xT))/ (xT^3) Sqrt[-1+E5t^2] x),{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,0.999999999999999,1},{E5t,1,(3-st)/(2 Sqrt[st])},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->20,WorkingPrecision->MachinePrecision]


averagesecondmom2[xT_?NumericQ,ac_,points_]:=-NIntegrate[(Sqrt[-1+cosh^2] Sqrt[9-4 st])/st^3 E^(-((3 cosh xT)/Sqrt[st])) Integrand[st,E5t,x,xs,\[Phi]] 1/(cosh E5t+Sqrt[-1+cosh^2] Sqrt[-1+E5t^2] x),{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,0.999999999999999,1},{E5t,1,(3-st)/(2 Sqrt[st])},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},{cosh,1,1.00000001,1.0001,1.01,1.1,1.2,1.3,1.5,2.,3.,5.,9,10,20,100.,1000.,100000.},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->20,WorkingPrecision->MachinePrecision]


averagesecondmom1[8.001939197344099,15,10^7]+averagesecondmom2[8.001939197344099,15,10^8]


averagesecondmomFUL[xT_?NumericQ,ac_,points_]:=(*1/(1536 m^2 \[Pi]^4)*Sqrt[3]3/m m m 1/(2(2Pi)^4) 9m^4 9m^2*)NIntegrate[(Sqrt[-1+cosh^2] Sqrt[9-4 st])/st^3 E^(-((3 cosh xT)/Sqrt[st])) Integrand[st,E5t,x,xs,\[Phi]](cosh E5t+Sqrt[-1+cosh^2] Sqrt[-1+E5t^2] x-1/(cosh E5t+Sqrt[-1+cosh^2] Sqrt[-1+E5t^2] x)),{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,0.999999999999999},{E5t,1,(3-st)/(2 Sqrt[st])},{x,-1,1},{xs,-1,1},{\[Phi],0,2Pi},{cosh,1,1.00000001,1.0001,1.01,1.1,1.2,1.3,1.5,2.,3.,5.,9,10,20,100.,1000.,100000.},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->20,WorkingPrecision->MachinePrecision]


averagesecondmomFUL[8.001939197344099,10,2*10^8]//Timing


1/(1536 m^2 \[Pi]^4)*Sqrt[3] 3/m m m 1/(2(2Pi)^4) 9m^4 9m^2


(* ::Subsection:: *)
(*Cannibal Reactions tabulated*)


TabulateP2Eaveraged[]:=Module[{
dataC0={},
dataC1={},
dataC2={},
xDMi=0.01,
xDMf=200,
Nn=120
},
Do[(*AppendTo[dataC0,{10^x,average2to3i[10^x]}];
Export["files/dataC0.dat",dataC0];*)
AppendTo[dataC1,{10^x,average2to3wP2Ei[10^x]}];
Export["files/dataC1V2.dat",dataC1];
(*AppendTo[dataC2,{10^x,averagesecondmomFUL[10^x,12,3*10^7]}];
Export["files/dataC2V2.dat",dataC2]*),
{x,Log10[xDMi]//N,Log10[xDMf]//N,(Log10[xDMf]-Log10[xDMi])/Nn//N}
]
]


TabulateP2Eaveraged[]


P2Eaveraged[xT_?NumericQ,points_?NumericQ,ac_?NumericQ]:=(*1/nphieq[xT,m]^39m^2mm^2/2*)(*factor*) NIntegrate[(E^(-((3 xT)/Sqrt[st]))  Sqrt[-(4/3)+3/st] Sqrt[-((-3+2 E5t Sqrt[st]+st)/st)] (E5t (2 (-1+E5t^2) st xA (Sqrt[st]+3 xT)+9 E^((3 xT)/Sqrt[st]) E5t Sqrt[(-1+E5t^2) st] xT^2 BesselK[2,(3 xT)/Sqrt[st]])))/(E5t Sqrt[9-6 E5t Sqrt[st]+st] xT^3) M[st,E5t,xA,xs,\[Phi]]^2 1/st^2,{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,1},{E5t,1,-( (-3+st)/(2 Sqrt[st]))},{xA,-1,0,1},{xs,-1,0,1},{\[Phi],0,2Pi},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->13,WorkingPrecision->MachinePrecision];


P2Eaverageds2[xT_?NumericQ,points_?NumericQ,ac_?NumericQ]:=(*1/nphieq[xT,m]^39m^2mm^2/2*)-(*factor*)  27NIntegrate[(*(  9 Sqrt[3]Sqrt[-1+E5t^2] Sqrt[-4+9/st] Sqrt[-((-3+2 E5t Sqrt[st]+st)/st)] )/Sqrt[9-6 E5t Sqrt[st]+st] (Sqrt[-1+z^2])/(E5t z+Sqrt[-1+E5t^2] xA Sqrt[-1+z^2])*)( E^(-((3 xT z)/Sqrt[st])) Sqrt[-1+E5t^2] Sqrt[3-(4 st)/3] Sqrt[3-2 E5t Sqrt[st]-st] Sqrt[-1+z^2])/(E5t st Sqrt[9-6 E5t Sqrt[st]+st] (z+(Sqrt[-1+E5t^2] xA Sqrt[-1+z^2])/E5t)) M[st,E5t,xA,xs,\[Phi]]^2 1/st^2,{st,10^-10,10^-5,10^-4,10^-2,10^-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99,0.999,1},{E5t,1,-( (-3+st)/(2 Sqrt[st]))},{xA,-1,0,1},{xs,-1,0,1},{\[Phi],0,2Pi},{z,1,1.01,1.1,1.2,1.3,2.,3.,5.,10.,20.,30.,100.},Method->{"QuasiMonteCarlo","MaxPoints"->points,"SymbolicProcessing"->0.},AccuracyGoal->ac,PrecisionGoal->13,WorkingPrecision->MachinePrecision];
