void MakeCDCConstraintsFile(bool perRing = true){

    // Some CDC geometry information
    int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};
    int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
    double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};
    double phi_ds[28] = { 0.000557611,0.0764693,0.0385138,0.0975182,-0.816345,-0.864077,-0.696401,-0.736506,0.656304,0.690227,0.57023,0.599326,0.0410675,0.0145592,0.00729358,0.0296972,0.43739,0.415211,0.38506,0.405461,-0.355973,-0.337391,-0.317012,-0.334703,0.0212654,0.0058214,0.005997358,0.0213175};
    // Four output files
    ofstream outFile_dxu;
    outFile_dxu.open("dxu_Constraints.txt");
    ofstream outFile_dxd;
    outFile_dxd.open("dxd_Constraints.txt");
    ofstream outFile_dyu;
    outFile_dyu.open("dyu_Constraints.txt");
    ofstream outFile_dyd;
    outFile_dyd.open("dyd_Constraints.txt");

    if (!perRing){
        outFile_dxu << "Contraint 0.0" << endl;
        outFile_dxd << "Contraint 0.0" << endl;
        outFile_dyu << "Contraint 0.0" << endl;
        outFile_dyd << "Contraint 0.0" << endl;
    }
    // No shrinking or stretching of CDC
    // Loop over the rings
    //ofstream outFile_dxu;
    //ofstream outFile_dxd;
    //ofstream outFile_dyu;
    //ofstream outFile_dyd;

    for (unsigned int iRing = 1; iRing <= 28; iRing++){
        //Get angular spacing between straws
        double dPhi = 2*TMath::Pi() / Nstraws[iRing-1];
        if (perRing){
            outFile_dxu << "Contraint 0.0" << endl;
            outFile_dxd << "Contraint 0.0" << endl;
            outFile_dyu << "Contraint 0.0" << endl;
            outFile_dyd << "Contraint 0.0" << endl;
        }
        //Loop over straws
        for (unsigned int iStraw = 1; iStraw <= Nstraws[iRing-1]; iStraw++){
            int index = straw_offset[iRing]+iStraw;
            int dxu_index = 1000 + (index-1)*4 +1;
            int dyu_index = 1000 + (index-1)*4 +2;
            int dxd_index = 1000 + (index-1)*4 +3;
            int dyd_index = 1000 + (index-1)*4 +4; 

            double ConstraintXU = TMath::Cos(phi[iRing-1] + dPhi*(iStraw-1));
            double ConstraintXD = TMath::Cos(phi_ds[iRing-1] + dPhi*(iStraw-1));
            double ConstraintYU = TMath::Sin(phi[iRing-1] + dPhi*(iStraw-1));
            double ConstraintYD = TMath::Sin(phi_ds[iRing-1] + dPhi*(iStraw-1));

            outFile_dxu << dxu_index << " " << ConstraintXU << endl;
            outFile_dxd << dxd_index << " " << ConstraintXD << endl;
            outFile_dyu << dyu_index << " " << ConstraintYU << endl;
            outFile_dyd << dyd_index << " " << ConstraintYD << endl;

        }
    }

    if (!perRing){
        outFile_dxu << "Contraint 0.0" << endl;
        outFile_dxd << "Contraint 0.0" << endl;
        outFile_dyu << "Contraint 0.0" << endl;
        outFile_dyd << "Contraint 0.0" << endl;
    }
    // No shrinking or stretching of CDC
    // Loop over the rings
    //ofstream outFile_dxu;
    //ofstream outFile_dxd;
    //ofstream outFile_dyu;
    //ofstream outFile_dyd;

    for (unsigned int iRing = 1; iRing <= 28; iRing++){
        //Get angular spacing between straws
        double dPhi = 2*TMath::Pi() / Nstraws[iRing-1];
        if (perRing){
            outFile_dxu << "Contraint 0.0" << endl;
            outFile_dxd << "Contraint 0.0" << endl;
            outFile_dyu << "Contraint 0.0" << endl;
            outFile_dyd << "Contraint 0.0" << endl;
        }
        //Loop over straws
        for (unsigned int iStraw = 1; iStraw <= Nstraws[iRing-1]; iStraw++){
            int index = straw_offset[iRing]+iStraw;
            int dxu_index = 1000 + (index-1)*4 +1;
            int dyu_index = 1000 + (index-1)*4 +2;
            int dxd_index = 1000 + (index-1)*4 +3;
            int dyd_index = 1000 + (index-1)*4 +4;

            double ConstraintXU = TMath::Sin(phi[iRing-1] + dPhi*(iStraw-1));
            double ConstraintXD = TMath::Sin(phi_ds[iRing-1] + dPhi*(iStraw-1));
            double ConstraintYU = TMath::Cos(phi[iRing-1] + dPhi*(iStraw-1));
            double ConstraintYD = TMath::Cos(phi_ds[iRing-1] + dPhi*(iStraw-1));

            outFile_dxu << dxu_index << " " << ConstraintXU << endl;
            outFile_dxd << dxd_index << " " << ConstraintXD << endl;
            outFile_dyu << dyu_index << " " << ConstraintYU << endl;
            outFile_dyd << dyd_index << " " << ConstraintYD << endl;
        }
    }
    /*
    // Suppress next order of harmonics
    if (!perRing){
    outFile_dxu << "Contraint 0.0" << endl;
    outFile_dxd << "Contraint 0.0" << endl;
    outFile_dyu << "Contraint 0.0" << endl;
    outFile_dyd << "Contraint 0.0" << endl;
    }
    for (unsigned int iRing = 1; iRing <= 28; iRing++){
//Get angular spacing between straws
double dPhi = 2*TMath::Pi() / Nstraws[iRing-1];
if (perRing){
outFile_dxu << "Contraint 0.0" << endl;
outFile_dxd << "Contraint 0.0" << endl;
outFile_dyu << "Contraint 0.0" << endl;
outFile_dyd << "Contraint 0.0" << endl;
}
//Loop over straws
for (unsigned int iStraw = 1; iStraw <= Nstraws[iRing-1]; iStraw++){
int index = straw_offset[iRing]+iStraw;
int dxu_index = (index-1)*4 +1;
int dxd_index = (index-1)*4 +2;
int dyu_index = (index-1)*4 +3;
int dyd_index = (index-1)*4 +4;

double ConstraintXU = TMath::Cos(phi[iRing-1] + dPhi*(iStraw-1)*2);
double ConstraintXD = TMath::Cos(phi_ds[iRing-1] + dPhi*(iStraw-1)*2);
double ConstraintYU = TMath::Sin(phi[iRing-1] + dPhi*(iStraw-1)*2);
double ConstraintYD = TMath::Sin(phi_ds[iRing-1] + dPhi*(iStraw-1)*2);

outFile_dxu << dxu_index << " " << ConstraintXU << endl;
outFile_dxd << dxd_index << " " << ConstraintXD << endl;
outFile_dyu << dyu_index << " " << ConstraintYU << endl;
outFile_dyd << dyd_index << " " << ConstraintYD << endl;

}
}
*/
outFile_dxu << "Contraint 0.0" << endl;
outFile_dxd << "Contraint 0.0" << endl;
outFile_dyu << "Contraint 0.0" << endl;
outFile_dyd << "Contraint 0.0" << endl;

//No global shifts
for (unsigned int iRing = 1; iRing <= 28; iRing++){
    //Get angular spacing between straws
    double dPhi = 2*TMath::Pi() / Nstraws[iRing-1];
    //Loop over straws
    for (unsigned int iStraw = 1; iStraw <= Nstraws[iRing-1]; iStraw++){
        int index = straw_offset[iRing]+iStraw;
        int dxu_index = 1000 + (index-1)*4 +1;
        int dyu_index = 1000 + (index-1)*4 +2;
        int dxd_index = 1000 + (index-1)*4 +3;
        int dyd_index = 1000 + (index-1)*4 +4;

        double ConstraintX = 1;
        double ConstraintY = 1;

        outFile_dxu << dxu_index << " " << ConstraintX << endl;
        outFile_dxd << dxd_index << " " << ConstraintX << endl;
        outFile_dyu << dyu_index << " " << ConstraintY << endl;
        outFile_dyd << dyd_index << " " << ConstraintY << endl;

    }
}

outFile_dxu.close();
outFile_dxd.close();
outFile_dyu.close();
outFile_dyd.close();
return;
}
