{
  // draws a single track as found in fdcData.txt

  Double_t PI = 3.1415962;

  // create a canvas
  TCanvas *universe = new TCanvas( "universe", "track", 100, 100, 700, 500 );

  // define the geometry
  gSystem->Load("libGeom"); // load geometry routines
  new TGeoManager("world", "the simplest geometry"); // define global
						     // geometry
						     // manager
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0,0); // define
							  // vacuum
							  // material
  TGeoMedium *med = new TGeoMedium("Vacuum", 1, mat); // define medium
						      // made of
						      // vacuum
  // define a box, destined to be the top volume
  TGeoVolume *top = gGeoManager->MakeBox("Top", med, 100.0, 100.0, 500.0);

  // define sphere shape
  TGeoSphere *sphere = new TGeoSphere(0.0, 0.2, 0.0, 180.0, 0.0, 360.0);
  // define the volume bubble as a sphere made of vacuum
  TGeoVolume *sphereVac = new TGeoVolume("bubble", sphere, med);
  sphereVac->SetLineColor(4);
  TGeoTranslation *trans[100]; // define an array of translations for hits
  Double_t x, y, z, theta, phi, dist; // to receive data from file
  Double_t thetadeg, phideg; // to hold angles in degrees
  Int_t index; // ditto
  ifstream inputFDC("fdcData.txt"); // open the file
  Double_t xmin = 1.e10; // initialize the bounding parameters
  Double_t xmax = -1.e10;
  Double_t ymin = 1.e10;
  Double_t ymax = -1.e10;
  Double_t zmin = 1.e10;
  Double_t zmax = -1.e10;
  Int_t nhitsFDC = 0; // zero hit counter
  while (inputFDC >> index >> x >> y >> z) { // while data can be read
					  // from the file
    cout << index << ' ' << x << ' ' << y << ' ' << z << endl;
    if (x < xmin) xmin = x; // test for new bounding parameters
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
    // store the translation for this hit
    trans[nhitsFDC] = new TGeoTranslation(x, y, z);
    top->AddNode(sphereVac, nhitsFDC, trans[nhitsFDC]); // add the ball to
						  // the world volume
    nhitsFDC++; // increment hit counter
  }
  inputFDC.close();
  /*
  cout << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax
       << ' ' << zmin << ' ' << zmax << endl;
  */

  TGeoTranslation transCDC; // translations for hits
  TGeoRotation rotCDC; // rotation
  TGeoCombiTrans *combiCDC[100]; // an array of combo tranlations and rotations
  ifstream inputCDC("cdcData.txt"); // open the file
  Int_t nhitsCDC = 0; // zero hit counter
  while (inputCDC >> index >> x >> y >> z >> theta >> phi >> dist) { // while data can be read
					  // from the file
    if (x < xmin) xmin = x; // test for new bounding parameters
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
    // store the translation for this hit
    transCDC.SetTranslation(x, y, z);
    phideg = phi*180.0/PI + 90;
    thetadeg = theta*180.0/PI;
    rotCDC.SetAngles(phideg, thetadeg, 0.0);
    cout << index << ' ' << x << ' ' << y << ' ' << z << ' ' << thetadeg << ' ' << phideg << endl;
    // define straw shape as a cylinder
    TGeoTube *straw = new TGeoTube(0.0, dist, 1.0);
    // define the straw as a cylinder made of vacuum
    TGeoVolume *strawVac = new TGeoVolume("straw", straw, med);
    if (abs(thetadeg) < 0.1) { // axial
      strawVac->SetLineColor(6);
    } else { // stereo
      strawVac->SetLineColor(7);
    }
    combiCDC[nhitsCDC] = new TGeoCombiTrans(transCDC, rotCDC);
    top->AddNode(strawVac, nhitsCDC, combiCDC[nhitsCDC]); // add the straw to
						  // the world volume
    nhitsCDC++; // increment hit counter
  }
  inputCDC.close();
  /*
  cout << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax
       << ' ' << zmin << ' ' << zmax << endl;
  */

  Double_t rmin[3], rmax[3]; // store bounds into an array
  rmin[0] = xmin - 1.0;
  rmax[0] = xmax + 1.0;
  rmin[1] = ymin - 1.0;
  rmax[1] = ymax + 1.0;
  rmin[2] = zmin - 1.0;
  rmax[2] = zmax + 1.0;

  TVirtualGeoTrack* track[2];
  Int_t track_index, track_id;
  for (track_id = 0; track_id < 2; track_id++) {
    track_index = gGeoManager->AddTrack(track_id,0,0); // create a track
    track[track_id] = gGeoManager->GetTrack(track_index); // get its pointer
  }
  //read in the track
  ifstream intraj("traj.txt"); // open the trajectory file
  Int_t npoints = 0; // zero point counter
  Int_t id; // track id (starting or final)
  while (intraj >> track_id >> index >> x >> y >> z) { // while data can be read
					  // from the file
    /*
    cout << track_id << ' ' << index << ' ' << x << ' ' << y << ' ' << z << endl;
    */
    // store the point on the track
    track[track_id]->AddPoint(x, y, z, (Double_t)(index)); // use index as a
						 // surrogate for time
    npoints++; // increment point counter
  }
  intraj.close();

  // define virtual track bounding box for visual effect only
  TGeoBBox *trackBox = new TGeoBBox("trackBox", 0.5*(xmax - xmin), 0.5*(ymax - ymin), 0.5*(zmax - zmin));
  TGeoVolume *trackVol = new TGeoVolume("track_volume", trackBox, med);
  TGeoTranslation *trackTrans = new TGeoTranslation((xmax + xmin)*0.5, (ymax + ymin)*0.5, (zmax + zmin)*0.5);
  top->AddNode(trackVol, 0, trackTrans);

  gGeoManager->SetTopVolume(top); // make top the top volume
  gGeoManager->CloseGeometry(); // we are finished with creating the geometry
  top->SetLineColor(kMagenta); // set color of the universe
  gGeoManager->SetTopVisible(); // allow top volume to be seen

  top->Draw(); // draw it
  track[0]->SetLineColor(1);
  track[1]->SetLineColor(2);
  track[0]->Draw();
  track[1]->Draw();

  // draw some axes
  TAxis3D *axes = new TAxis3D();
  axes->Draw();

  // zoom in on the track
  TView *view = TView::CreateView(1, rmin, rmax);

  universe->Update(); 
}
