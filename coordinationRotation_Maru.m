%NVec: unit vector representing the microcolumn axis or apical dendrites
ClmAxisVec=NVec; 

ClmAxisVec(3)=sqrt(1-(ClmAxisVec(1)^2+ClmAxisVec(2)^2));
[tgtrng,tgtang] = rangeangle(ClmAxisVec');

TempAngle=[tgtang(1),tgtang(2)-90];

Theta2=deg2rad(-1*TempAngle(1));
Phi2=deg2rad(-1*TempAngle(2));
ClmPrjMtx1=[cos(Theta2) -sin(Theta2);sin(Theta2) cos(Theta2)];
ClmPrjMtx2=[cos(Phi2) -sin(Phi2);sin(Phi2) cos(Phi2)];

TempVec=ClmAxisVec;

Temp4=ClmPrjMtx1*[TempVec(1);TempVec(2)];
TempVec(1)=Temp4(1);
TempVec(2)=Temp4(2);
Temp3=ClmPrjMtx2*[TempVec(1);TempVec(3)];
TempVec(1)=Temp3(1);
TempVec(3)=Temp3(2);
[temprng,tempang] = rangeangle(TempVec');

disp('ClmAxisVecの回転後の行列')
TempVec

%%%

TempPosUM1=Pos; %Pos: 3D coordinates of GCaMP6 cells


tTempPosUM1=TempPosUM1;
h = waitbar(0,'--座標回転中1--');
for iCell=1:size(TempPosUM1,1)
    Temp1=ClmPrjMtx1*TempPosUM1(iCell,1:2)';
    tTempPosUM1(iCell,1)=Temp1(1);
    tTempPosUM1(iCell,2)=Temp1(2);
    Temp2=ClmPrjMtx2*[tTempPosUM1(iCell,1);TempPosUM1(iCell,3)];
    tTempPosUM1(iCell,1)=Temp2(1);
    tTempPosUM1(iCell,3)=Temp2(2);

    waitbar(iCell/size(tTempPosUM1,1));
end
close(h)