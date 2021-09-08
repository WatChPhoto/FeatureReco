void dhist(){
  TTree t_photogram("t_photogram", "Photogrammetry Results");
  t_photogram.ReadFile("SK_ring_cameras.txt");

  Char_t          CameraID[4];
  Double_t        CameraPosition[3];
  Double_t        CameraOrientation[3][3];
  t_photogram.SetBranchAddress("CameraID", CameraID);
  t_photogram.SetBranchAddress("CameraPosition", CameraPosition);
  t_photogram.SetBranchAddress("CameraOrientation", CameraOrientation);
  
  int n=t_photogram.GetEntries();
  
  std::cout<<n<<std::endl;
  for(int i=0;i<n;i++){
    t_photogram.GetEntry(i);
    std::cout<<CameraID<<" "<<CameraOrientation[0][0]<<" "<<CameraOrientation[0][1]<<" "<<CameraOrientation[0][2]<<" "<<CameraOrientation[1][0]<<" "<<CameraOrientation[1][1]<<" "<<CameraOrientation[1][2]<<" "<<CameraOrientation[2][0]<<" "<<CameraOrientation[2][1]<<" "<<CameraOrientation[2][2]<<std::endl;
  }
}
