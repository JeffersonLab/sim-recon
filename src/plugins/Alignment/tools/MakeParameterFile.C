void MakeParameterFile(){

   double scale = 1.0;
   double translationPresigma = 0.001;
   double rotationPresigma = 0.00001;
   double gainPresigma = 0.01;

   ofstream outfile;
   outfile.open("Parameters.txt");

   bool fixGains = true;

   outfile << "Parameter" << endl;

   // FDC
   for (unsigned int i=1; i<=24; i++){
      
      int indexOffset = 100000 + i*1000;
      // dXw, dYw
      outfile << indexOffset + 1 << " 0.0 " << translationPresigma << endl;
      outfile << indexOffset + 100 << " 0.0 " << translationPresigma << endl;

      // dU, dV
      outfile << indexOffset + 101 << " 0.0 " << translationPresigma << endl;
      outfile << indexOffset + 102 << " 0.0 " << translationPresigma << endl;

      // dPhiU, dPhiV
      outfile << indexOffset + 103 << " 0.0 " << rotationPresigma << endl;
      outfile << indexOffset + 104 << " 0.0 " << rotationPresigma << endl;

      // dPhiWz, dPhiWx
      outfile << indexOffset + 2 << " 0.0 0.001" << endl;
      outfile << indexOffset + 3 << " 0.0 -1.0 " << endl; // Fix for now

      //Strip Pitch
      for (unsigned int j = 200; j <= 209; j++){
         outfile << indexOffset + j << " 0.0 " << translationPresigma << endl;
      }

      //Strip Gain
      for (unsigned int j = 301; j <= 516; j++){
         if (fixGains) outfile << indexOffset + j << " 0.0 -1.0" << endl;
         else outfile << indexOffset + j << " 0.0 " << gainPresigma << endl;
      }

      for (unsigned int j = 601; j <= 816; j++){
         if (fixGains) outfile << indexOffset + j << " 0.0 -1.0" << endl;
         else outfile << indexOffset + j << " 0.0 " << gainPresigma << endl;
      }
   }

   outfile.close();
}


