(* ::Package:: *)

BeginPackage@"Functions`";
If[Position[$Path,NotebookDirectory[]]==={},AppendTo[$Path,NotebookDirectory[]]];

matrixZRotation::usage="Gives rotation matrix around Z axis. Accepts rotatoin angle.";
matrixXRotation::usage="Gives rotation matrix around X axis. Accepts rotatoin angle.";
matrixTranslation::usage="Gives translation matrix. Accepts numbers: x, y, z.";
inertiaTensor::usage="Calculates inertia tensor.";
dynamicEquations::usage="Calculate dynamics with Newton-Euler algorithm.";

Begin@"`Private`";

matrixZRotation[theta_]:={
	{Cos@theta,-Sin@theta,0,0},
	{Sin@theta,Cos@theta,0,0},
	{0,0,1,0},
	{0,0,0,1}};
matrixXRotation[theta_]:={
	{1, 0, 0, 0},
	{0,Cos@theta,-Sin@theta,0},
	{0,Sin@theta,Cos@theta,0},
	{0,0,0,1}};
matrixTranslation[x_,y_,z_]:={
	{1,0,0,x},
	{0,1,0,y},
	{0,0,1,z},
	{0,0,0,1}};
inertiaTensor[xFrom_,xTo_,yFrom_,yTo_,zFrom_,zTo_,p_]:=Module[{x,y,z,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,integral},
	integral[f_]:=Integrate[f p,{z,zFrom,zTo},{y,yFrom,yTo},{x,xFrom,xTo}];
	Ixx=integral[y^2+z^2];
	Iyy=integral[x^2+z^2];
	Izz=integral[x^2+y^2];
	Ixy=integral[x y];
	Ixz=integral[x z];
	Iyz=integral[y z];
	{{Ixx,-Ixy,-Ixz},{-Ixy,Iyy,-Ixz},{-Ixz,-Iyz,Izz}}//Simplify
];
dynamicEquations[Ts_,Is_,ms_,Ps_,fn_,nn_,w0_,dw0_,dv0_]:=Module[{nodes,ws,dws,dvs,dvcs,Fs,Ns,fs,ns,ts,getR,R,getP,P,thetaD,thetaDD,thetaDZ,thetaDDZ},
	getR[T_]:=Take[T,{1,3},{1,3}];
	getP[T_]:=Take[T,{1,3},{4,4}]//Flatten;
	nodes=(Length@Ts)+1;
	ws={};
	dws={};
	dvs={};
	dvcs={};
	Fs={};
	Ns={};
	fs={};
	ns={};
	ts={};
	Do[
		ws~AppendTo~{0,0,0};
		dws~AppendTo~{0,0,0};
		dvs~AppendTo~{0,0,0};
		dvcs~AppendTo~{0,0,0};
		Fs~AppendTo~{0,0,0};
		Ns~AppendTo~{0,0,0};
		fs~AppendTo~{0,0,0};
		ns~AppendTo~{0,0,0};
		ts~AppendTo~{0,0,0};
	,{i,1,nodes+1}];
	fs[[nodes]]=fn;
	ns[[nodes]]=nn;
	ws[[1]]=w0;
	dws[[1]]=dw0;
	dvs[[1]]=dv0;
	Do[
		R=getR@Ts[[i]]//Transpose;
		P=getP@Ts[[i]];
		thetaDZ={0,0,Symbol@StringJoin["theta",ToString@i,"D"]};
		thetaDDZ={0,0,Symbol@StringJoin["theta",ToString@i,"DD"]};
		ws[[i+1]]=R.ws[[i]]+thetaDZ;
		dws[[i+1]]=R.dws[[i]]+Cross[R.ws[[i]],thetaDZ]+thetaDDZ;
		dvs[[i+1]]=R.(Cross[dws[[i]],P]+Cross[ws[[i]],Cross[ws[[i]],P]]+dvs[[i]]);
		dvcs[[i+1]]=Cross[dws[[i+1]],Ps[[i]]]+Cross[ws[[i+1]],Cross[ws[[i+1]],Ps[[i]]]]+dvs[[i+1]];
		Fs[[i+1]]=ms[[i]]*dvcs[[i+1]];
		Ns[[i+1]]=Is[[i]].dws[[i+1]]+Cross[ws[[i+1]],Is[[i]].ws[[i+1]]];
	,{i,1,nodes-1}];
	Do[
		R=getR@Ts[[i-1]];
		P=getP@Ts[[i-1]];
		fs[[i]]=R.fs[[i+1]]+Fs[[i]];
		ns[[i]]=Ns[[i]]+R.ns[[i+1]]+Cross[Ps[[i-1]],Fs[[i]]]+Cross[P,R.fs[[i+1]]];
		ts[[i]]=ns[[i]].{0,0,1};
	,{i,nodes,2,-1}];
	ts
]; 

End[];
EndPackage[];
