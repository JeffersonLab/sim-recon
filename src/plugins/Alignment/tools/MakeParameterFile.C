void MakeParameterFile(){

   // Choose which parameters you would like to fix
   bool fixCDCWires = true;
   bool fixCDCt0 = true;
   bool fixCDCGlobalX = true;
   bool fixCDCGlobalY = true;
   bool fixCDCGlobalZ = true;
   bool fixCDCGlobalPhiX = true;
   bool fixCDCGlobalPhiY = true;
   bool fixCDCGlobalPhiZ = true;
   bool fixFDCt0 = true;
   bool fixFDCGains = true;
   bool fixFDCCathodeOffsets = true;
   bool fixFDCCathodeAngles = true;
   bool fixFDCCellOffsetsWires = true;
   bool fixFDCCellOffsetsCathodes = true;
   bool fixFDCWireRotationX = true;
   bool fixFDCWireRotationY = true;
   bool fixFDCWireRotationZ = true;
   bool fixFDCZ = true;
   bool fixFDCPitch = true;
   bool fixFDCGap = true;

   double scale = 1.0;
   double translationPresigma = 0.0005;
   double rotationPresigma = 0.0001;
   double gainPresigma = 0.01;
   double pitchPresigma = 0.001;
   double t0Presigma=10.0;

   ofstream outfile;
   outfile.open("Parameters.txt");

   outfile << "Parameter" << endl;

   // CDC
   if(fixCDCGlobalX) outfile << "1 0.0 -1" << endl;
   else outfile << "1 0.0 " << translationPresigma << endl;

   if(fixCDCGlobalY) outfile << "2 0.0 -1" << endl;
   else outfile << "2 0.0 " << translationPresigma << endl;

   if(fixCDCGlobalZ) outfile << "3 0.0 -1" << endl;
   else outfile << "3 0.0 " << translationPresigma << endl;

   if(fixCDCGlobalPhiX) outfile << "4 0.0 -1" << endl;
   else outfile << "4 0.0 " << rotationPresigma << endl;

   if(fixCDCGlobalPhiY) outfile << "5 0.0 -1" << endl;
   else outfile << "5 0.0 " << rotationPresigma << endl;

   if(fixCDCGlobalPhiZ) outfile << "6 0.0 -1" << endl;
   else outfile << "6 0.0 " << rotationPresigma << endl;

   for (unsigned int i=16001; i <= 19522; i++){
      if (fixCDCt0) outfile << i << " 0.0 -1.0" << endl;
      else outfile << i << " 0.0 " << t0Presigma << endl;
   }

   for (unsigned int i=1001; i <= 15088; i++){
      if (fixCDCWires) outfile << i << " 0.0 -1.0" << endl;
      else outfile << i << " 0.0 " << translationPresigma << endl;
   }

   // Possibility of varying presigma around the detector
   /*
      int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};
      int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
      double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
      double phiOffset[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};

   // Need to use smaller presigma for values that are not very aligned for
   unsigned int ring = 1;
   for (unsigned int i=1001; i <= 15088; i++){
   unsigned int strawNumber = (i-1001)/4 +1;
   strawNumber-=straw_offset[ring];
   double phi = phiOffset[ring-1] + (strawNumber -1) * (2*TMath::Pi()/Nstraws[ring-1]);
   // sin for x cos for y
   double xTranslationPresigma = translationPresigma - (0.009) * cos(phi) * cos(phi);
   double yTranslationPresigma = translationPresigma - (0.009) * sin(phi) * sin(phi);
   // x offsets
   if (i%4 == 1 || i%4 == 3){
   outfile << i << " 0.0 " << xTranslationPresigma << endl;
   }
   //y offsets
   else{
   outfile << i << " 0.0 " << yTranslationPresigma << endl;
   }
   if (i%4 == 0 && strawNumber==Nstraws[ring-1]) ring++; 
   }
   */

   // FDC
   for (unsigned int i=1; i<=24; i++){

      int indexOffset = 100000 + i*1000;
      if(fixFDCZ){
         outfile<< indexOffset + 5 << " 0.0 -1.0" << endl;
      }
      else {
         outfile<< indexOffset + 5 << " 0.0 " << translationPresigma << endl;
      }
      if (fixFDCCellOffsetsWires){
         outfile << indexOffset + 1 << " 0.0 -1.0" << endl;
      }
      else{
         outfile << indexOffset + 1 << " 0.0 " << translationPresigma << endl;
      }
      if (fixFDCCellOffsetsCathodes){
         outfile << indexOffset + 100 << " 0.0 -1.0" << endl;
      }
      else{
         outfile << indexOffset + 100 << " 0.0 " << translationPresigma << endl;
      }

      if (fixFDCWireRotationX) outfile << indexOffset + 2 << " 0.0 -1.0" << endl;
      else outfile << indexOffset + 2 << " 0.0 " << rotationPresigma << endl;

      if (fixFDCWireRotationY) outfile << indexOffset + 3 << " 0.0 -1.0" << endl;
      else outfile << indexOffset + 3 << " 0.0 " << rotationPresigma << endl;

      if (fixFDCWireRotationZ) outfile << indexOffset + 4 << " 0.0 -1.0" << endl;
      else outfile << indexOffset + 4 << " 0.0 " << rotationPresigma << endl;

      if (fixFDCCathodeOffsets){
         outfile << indexOffset + 101 << " 0.0 -1.0" << endl;
         outfile << indexOffset + 102 << " 0.0 -1.0" << endl;
      }
      else{
         outfile << indexOffset + 101 << " 0.0 " << translationPresigma << endl;
         outfile << indexOffset + 102 << " 0.0 " << translationPresigma << endl;
      }

      if (fixFDCCathodeAngles){
         // dPhiU, dPhiV
         outfile << indexOffset + 103 << " 0.0 -1.0" << endl;
         outfile << indexOffset + 104 << " 0.0 -1.0" << endl;
      }
      else{
         // dPhiU, dPhiV
         outfile << indexOffset + 103 << " 0.0 " << rotationPresigma << endl;
         outfile << indexOffset + 104 << " 0.0 " << rotationPresigma << endl;
      }

      //Strip Pitch
      if (fixFDCPitch){
         outfile <<  indexOffset + 200 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 202 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 204 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 205 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 207 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 209 << " 0.0 -1.0" << endl;
      }
      else{
         outfile <<  indexOffset + 200 << " 0.0 " << pitchPresigma << endl;
         outfile <<  indexOffset + 202 << " 0.0 " << pitchPresigma << endl;
         outfile <<  indexOffset + 204 << " 0.0 " << pitchPresigma << endl;
         outfile <<  indexOffset + 205 << " 0.0 " << pitchPresigma << endl;
         outfile <<  indexOffset + 207 << " 0.0 " << pitchPresigma << endl;
         outfile <<  indexOffset + 209 << " 0.0 " << pitchPresigma << endl;
      }
      //Strip gap

      if (fixFDCGap){
         outfile <<  indexOffset + 201 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 203 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 206 << " 0.0 -1.0" << endl;
         outfile <<  indexOffset + 208 << " 0.0 -1.0" << endl;
      }
      else{
         outfile <<  indexOffset + 201 << " 0.0 " << translationPresigma << endl;
         outfile <<  indexOffset + 203 << " 0.0 " << translationPresigma << endl;
         outfile <<  indexOffset + 206 << " 0.0 " << translationPresigma << endl;
         outfile <<  indexOffset + 208 << " 0.0 " << translationPresigma << endl;
      }

      //Strip Gain
      for (unsigned int j = 301; j <= 517; j++){
         if (fixFDCGains) outfile << indexOffset + j << " 0.0 -1.0" << endl;
         else outfile << indexOffset + j << " 0.0 " << gainPresigma << endl;
      }

      for (unsigned int j = 601; j <= 817; j++){
         if (fixFDCGains) outfile << indexOffset + j << " 0.0 -1.0" << endl;
         else outfile << indexOffset + j << " 0.0 " << gainPresigma << endl;
      }

      //t0
      for (unsigned int j = 901; j <= 996; j++){
         if (fixFDCt0) outfile << indexOffset + j << " 0.0 -1.0" << endl;
         else outfile << indexOffset + j << " 0.0 " << t0Presigma << endl;
      }

   }

   outfile.close();
}


