(* ::Package:: *)

(* ::Input:: *)
g=MatrixPlot[Table[i + j*3, {i, 10}, {j, 10}], ColorFunction -> "TemperatureMap", ImageSize -> Large, ColorFunctionScaling -> True, LabelStyle -> Directive[Blue, FontFamily -> "Latin Modern Roman", FontSize -> 24]]
Export["testing_1.png", g, "PNG", ImageResolution->200];
