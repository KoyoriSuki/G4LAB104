


void readResults()
{
	const char* files[] = { "sh0.root", "sh50.root", "sh100.root", "sh120.root", "sh140.root" };
	const int nFiles = 5;

	// Create Canvas
	TCanvas* c1 = new TCanvas("c1", "Data Comparison", 800, 600);

	// Create Legend
	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

	for (int i = 0; i < nFiles; ++i)
	{
		TFile* file = TFile::Open(files[i]);
		TTree* tree = (TTree*)file->Get("TreeOfData");

		vector<double>* zenith_Angle = 0;
		vector<double>* initial_energy = 0;
		tree->SetBranchAddress("Direction_Zenith_Angle", &zenith_Angle);
		tree->SetBranchAddress("Initial_Energy", &initial_energy);
		Long64_t n = tree->GetEntries();
		for (Long64_t i = 0; i < n; i++) {
			tree->GetEntry(i);
			ULong_t nsig = zenith_Angle->size(); // MUST be done for every tree entry
			double* sig_zenith_Angle = zenith_Angle->data(); // MUST be done for every tree entry
			double* sig_initial_energy = initial_energy->data(); // MUST be done for every tree entry
			for (ULong_t j = 0; j < nsig; j++) {
				std::cout << i << " : " << j << " : " << sig_zenith_Angle[j] << std::endl;
				std::cout << i << " : " << j << " : " << sig_initial_energy[j] << std::endl;
			}
		}
		tree->ResetBranchAddresses(); // disconnect from local variables
		delete zenith_Angle; // cleanup
		delete initial_energy;
		break;
	}

}