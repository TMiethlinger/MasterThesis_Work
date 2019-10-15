directorypath="/home/k3501/k354524/MasterThesis/FluidizedBed_RP_toTMiethlinger/DEM/post_60k_1/elements/";
numbersInString[s_]:=ToExpression@StringCases[s,NumberString];
filenames=SortBy[FileNames["*",directorypath],numbersInString][[2;;]];


MyPlot[data_]:=MatrixPlot[data,PlotRangePadding->0,PlotTheme->"Scientific",LabelStyle->Directive[Black,FontColor->Black,24,FontFamily->"CMU Serif",Italic],
FrameStyle->Directive[Black,FontColor->Black,20,FontFamily->"CMU Serif",Plain],FrameLabel->{"n","m"},FrameTicks->{{Automatic,None},{Automatic,None}},
ImageSize->{1000,1000},DataReversed->{True,False},ColorFunction->(RGBColor[1-#,1-#,1-#]&),ColorFunctionScaling->False];
Iterator:={x@#,xmin@#+dx2@#,xmax@#-dx2@#,dx@#}&;
ShiftToGrid[\[Zeta]_,{dx_,dx2_}]:=Floor[\[Zeta],dx]+dx2;
Distance[countgridi_,countgridj_]:=Total[(N[countgridi]-N[countgridj])^2];

dt=10;
tmin=1;
tmax=Length@filenames;

xmin[1]=-0.04;
xmax[1]=0.04;
xmin[2]=-0.0075;
xmax[2]=0.0075;
xmin[3]=0.0;
xmax[3]=0.25;

dx[1]=0.01;
dx[2]=0.0075/2;
dx[3]=0.0125/2;

Map[(dx2[#]=dx[#]/2)&,Range[1,3]];
Map[(d[#]={dx[#],dx2[#]})&,Range[1,3]];

gridpoints=Flatten[Apply[Table,Join[{Array[x,3]},Map[Iterator,Range[3]]]],2];

t=AbsoluteTime[];

Do[

datai=Map[Internal`StringToDouble[#]&,Map[StringSplit[#," "]&,ReadList[filenames[[ti]],String][[10;;All]]][[All,2;;4]],{2}];
datagridi=MapIndexed[ShiftToGrid[#1,d[#2[[2]]]]&,datai,{2}];
countgridi=Map[Count[datagridi,#]&,gridpoints];

Do[

dataj=Map[Internal`StringToDouble[#]&,Map[StringSplit[#," "]&,ReadList[filenames[[tj]],String][[10;;All]]][[All,2;;4]],{2}];
datagirdj=MapIndexed[ShiftToGrid[#1,d[#2[[2]]]]&,dataj,{2}];
countgridj=Map[Count[datagirdj,#]&,gridpoints];

dm[ti,tj]=Distance[countgridi,countgridj];

,{tj,ti+dt,tmax,dt}];

,{ti,tmin,tmax-dt,dt}];

t=AbsoluteTime[]-t;

Do[
dm[ti,tj]=If[ti==tj,0,If[ti<tj,dm[ti,tj],dm[tj,ti]]];
,{ti,tmin,tmax,dt},{tj,tmin,tmax,dt}];


distancematrix=Table[dm[ti,tj], {ti, tmin, tmax, dt}, {tj, tmin, tmax, dt}];
(* max=Max@Flatten@distancematrix; *)
(* distancematrix=distancematrix/max; *)
Export["/home/k3501/k354524/MasterThesis/Programming/FieldRecurrences/rho/distancematrix_rho_"<>ToString[Length[datai]]<>"_"<>ToString[dt]<>".txt", distancematrix, "Table", "FieldSeparators" -> " "];
(* myplot=MyPlot@distancematrix; *)
(* Export["/home/k3501/k354524/MasterThesis/Programming/FieldRecurrences/rho/distancematrix_rho_1.png", myplot, "PNG"]; *)
Export["/home/k3501/k354524/MasterThesis/Programming/FieldRecurrences/rho/time_"<>ToString[Length[datai]]<>"_"<>ToString[dt]<>".txt", ToString[t]];
