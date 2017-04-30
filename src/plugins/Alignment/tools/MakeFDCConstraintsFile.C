void MakeFDCConstraintsFile(TString rootFile = "hd_root.root"){

   // Select constraints
   bool constrainZ = false;
   bool constrainWireCathodeAlignment = false;
   bool constraint0 = false;

   ofstream outfile;
   outfile.open("FDCConstraints.txt");

   if(constrainZ){
      for (unsigned int i=1; i<=4; i++){
         outfile << "Constraint 0.0" << endl;
         for (unsigned int j=1;j<=6;j++){
            int indexOffset = 100000 + ((i-1)*6+j)*1000;
            outfile << indexOffset+5 << " 1.0" << endl;
         }
      }
   }

   if(constraint0){
      //t0
      outfile <<"Constraint 0.0" << endl;
      for (unsigned int i=1; i<=24; i++){ 
         int indexOffset = 100000 + i*1000;
         for(unsigned int j=901; j<=996;j++){
            outfile << indexOffset + j << " 1.0 " << endl;
         }
      }
   }

   for (unsigned int i=1; i<=24; i++){
      int indexOffset = 100000 + i*1000;
      int histIndex = i*1000;

      double angleOffset = ((i-1)%6)*TMath::Pi()/3.;
      outfile << "Measurement 0.0 0.02" << endl;
      outfile << indexOffset+1 << " " << sin(angleOffset) << endl;
      outfile << indexOffset+100 << " " << sin(angleOffset+TMath::Pi()/2.) << endl;
      outfile << "Measurement 0.0 0.02" << endl;
      outfile << indexOffset+1 << " " << cos(angleOffset) << endl;
      outfile << indexOffset+100 << " " << cos(angleOffset+TMath::Pi()/2.) << endl;
   }

   if(constrainWireCathodeAlignment){
      // Here we need to existing values to write the constraint
      TFile *file = TFile::Open(rootFile);
      TProfile *hFDCConstants;
      file->GetObject("AlignmentConstants/FDCAlignmentConstants",hFDCConstants);

      for (unsigned int i=1; i<=24; i++){
         // If this is true the cathode pitch is constrained such that the average pitch is unchanged
         // but each of the sections can float. If this is false, the pitch is adjusted the same ammount for
         // each foil.
         bool averagePitch = true;

         int indexOffset = 100000 + i*1000;
         int histIndex = i*1000;
         double phiu = TMath::DegToRad()*75.0;
         double phiv = TMath::Pi()-phiu;

         phiu+=hFDCConstants->GetBinContent(histIndex+103);
         phiv+=hFDCConstants->GetBinContent(histIndex+104);

         double pu1 = hFDCConstants->GetBinContent(histIndex+200);
         double pu2 = hFDCConstants->GetBinContent(histIndex+202);
         double pu3 = hFDCConstants->GetBinContent(histIndex+204);

         double pv1 = hFDCConstants->GetBinContent(histIndex+205);
         double pv2 = hFDCConstants->GetBinContent(histIndex+207);
         double pv3 = hFDCConstants->GetBinContent(histIndex+209);

         double pu = (48.*(pu1+pu3)+96.*pu2)/192.; // Average strip pitch
         double pv = (48.*(pv1+pv3)+96.*pv2)/192.; // Average strip pitch

         double sinphiu = sin(phiu);
         double sinphiv = sin(phiv);
         double sinphiumphiv = sin(phiu-phiv);
         double cosphiumphiv = cos(phiu-phiv);

         if(averagePitch){
            //Avg = 0
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 200 << " 48.0 " << endl;
            outfile << indexOffset + 202 << " 96.0 " << endl;
            outfile << indexOffset + 204 << " 48.0 " << endl;
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 205 << " 48.0 " << endl;
            outfile << indexOffset + 207 << " 96.0 " << endl;
            outfile << indexOffset + 209 << " 48.0 " << endl;
         }
         else{
            //Equal
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 200 << " 1.0 " << endl;
            outfile << indexOffset + 202 << " -1.0 " << endl;
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 200 << " 1.0 " << endl;
            outfile << indexOffset + 204 << " -1.0 " << endl;
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 205 << " 1.0 " << endl;
            outfile << indexOffset + 207 << " -1.0 " << endl;
            outfile << "Constraint 0.0" << endl;
            outfile << indexOffset + 205 << " 1.0 " << endl;
            outfile << indexOffset + 209 << " -1.0 " << endl;
         }
         // Constrain Wire/cathode alignmnet
         outfile << "Constraint 0.0" << endl;
         //outfile << "Measurement 0.0 0.00001" << endl;
         outfile << indexOffset + 200 << " " <<  sinphiv  << endl;
         outfile << indexOffset + 205 << " " <<  sinphiu <<endl;
         outfile << indexOffset + 103 << " " << -1*(pv+pu*cosphiumphiv)*sinphiv/sinphiumphiv << endl;
         outfile << indexOffset + 104 << " " << (pu+pv*cosphiumphiv)*sinphiu/sinphiumphiv << endl;
         outfile << "Constraint 0.0" << endl;
         //outfile << "Measurement 0.0 0.00001" << endl;
         outfile << indexOffset + 200 << " " <<  sinphiv  << endl;
         outfile << indexOffset + 205 << " " <<  -sinphiu <<endl;
         outfile << indexOffset + 103 << " " << (pv-pu*cosphiumphiv)*sinphiv/sinphiumphiv << endl;
         outfile << indexOffset + 104 << " " << (pu-pv*cosphiumphiv)*sinphiu/sinphiumphiv << endl;
      }
   }

      outfile.close();
   }


