void MakeFDCConstraintsFile(){

   ofstream outfile;
   outfile.open("FDCConstraints.txt");

   // FDC
   for (unsigned int i=1; i<=24; i++){
      
      int indexOffset = 100000 + i*1000;

      // dU, dV
      outfile << "Constraint 0.0" << endl;
      outfile << indexOffset + 101 << " 1.0 " << endl;
      outfile << indexOffset + 102 << " -1.0 " << endl;


      // dPhiU, dPhiV
      outfile << "Constraint 0.0" << endl;
      outfile << indexOffset + 103 << " 1.0 " << endl;
      outfile << indexOffset + 104 << " 1.0 " << endl;
   }

   outfile.close();
}


