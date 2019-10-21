#!/apps/mathematica-12.0.0/bin/math -script

(*
Start example:
./distance_matrix.wls -relfolderpath "/master_thesis_work/results/field_distance_matrix/Liggghts/N_60000_U_1.333_1_100_12000_8_8_8" -filename "field_distance_matrix.txt" -nticks 4 -tmax 5
*)

filename = $CommandLine[[Position[$CommandLine, "-filename"][[-1, 1]] + 1]];
relfolderpath = $CommandLine[[Position[$CommandLine, "-relfolderpath"][[-1, 1]] + 1]];
filepathname = FileNameJoin[{$HomeDirectory, relfolderpath, filename}];
(* Print[filepathname]; *)
data = Import[filepathname, "Table"];
(* Print[$CommandLine] *)
(* Print[filepathname]; *)

nticks = ToExpression[$CommandLine[[Position[$CommandLine, "-nticks"][[-1, 1]] + 1]]];
tmax = N@ToExpression[$CommandLine[[Position[$CommandLine, "-tmax"][[-1, 1]] + 1]]];

n=Length@data;
min=Min[data];
max=Max[data];
datanorm=(data - min) / (max-min);

xl=Row[{DisplayForm[SubscriptBox[Style["t", FontFamily->"Latin Modern Roman", Italic, 24], Style["1", 16, FontFamily->"Latin Modern Roman", Plain]]], Style[" [s]", FontFamily->"Latin Modern Roman", Plain, 24]}];
yl=Row[{DisplayForm[SubscriptBox[Style["t", FontFamily->"Latin Modern Roman", Italic, 24], Style["2", 16, FontFamily->"Latin Modern Roman", Plain]]], Style[" [s]", FontFamily->"Latin Modern Roman", Plain, 24]}];
frameticks = Table[{i + 0.5, tmax*(i/n)}, {i, 0, n, n/nticks}];

g=Row[{
MatrixPlot[
datanorm, 
PlotRange->All, 
DataReversed->{True, False}, 
ColorFunction->"TemperatureMap", 
ColorFunctionScaling->False, 
PlotTheme->"Scientific", 
FrameStyle->Directive[Black, FontColor->Black, 18, FontFamily->"Latin Modern Roman", Plain], 
FrameTicks -> {{frameticks, None}, {frameticks, None}},
FrameLabel -> {yl, xl},
ImagePadding->{{70, 15}, {70, 10}}, 
ImageSize->{Automatic, 500}, 
PlotRangePadding->0, 
AspectRatio->1.0
], 
DensityPlot[y, {x, -1, 1}, {y, -1, 1}, 
PlotRange->{{0, 1}, {0, 1}}, 
PlotPoints->100, 
MaxRecursion->3, 
ColorFunction->"TemperatureMap", 
ColorFunctionScaling->False, 
FrameStyle->Directive[Black, FontColor->Black, 18, FontFamily->"Latin Modern Roman", Plain], 
FrameTicks->{{None, All}, {None, None}}, 
ImagePadding->{{3, 40}, {70, 10}}, 
ImageSize->{Automatic, 500}, 
PlotRangePadding->0, 
AspectRatio->16]
}];

outputfolderpath = FileNameJoin[{$HomeDirectory, StringReplace[relfolderpath, "results" -> "figures"]}];
(* Print[outputfolderpath]; *)
Quiet[CreateDirectory[outputfolderpath]];
outputfilename = StringReplace[filename, ".txt" -> ".png"];
filepathname = FileNameJoin[{outputfolderpath, outputfilename}];

Quiet[DeleteFile[filepathname]];
Export[filepathname, g, "PNG", ImageResolution -> 500];

(* For testing purposes *)
(* n=200; *)
(* data=Table[Cos[2 \[Pi] 4/n i] Sin[2 \[Pi] 3/n j], {i, n}, {j, n}]; *)
