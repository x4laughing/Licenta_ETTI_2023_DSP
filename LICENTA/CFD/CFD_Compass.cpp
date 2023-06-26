

{
	std::vector<double> CFD(const std::vector<double> &rawSample, const std::vector<double> &sampleNumber, double baseline, float fraction_cfd, int delay_cfd)
	{
		std::vector<double> cfd;
		std::vector<double> cfd_att;
		std::vector<double> cfd_delay;

		for (auto i = 0; i < delay_cfd; i++)
		{
			cfd_delay.push_back(0);
		}

		for (auto &sig_iter : rawSample)
		{
			cfd_att.push_back((sig_iter - baseline) * fraction_cfd);
			cfd_delay.push_back(-sig_iter + baseline);
		}

		for (int i = 0; i < sampleNumber.size(); i++)
		{
			cfd.push_back(cfd_att[i] + cfd_delay[i]);
		}

		return cfd;
	}
	TFile *f_positron = TFile::Open("/home/georgen24/Documents/LICENTA/coduri/Root_files/run_same_settings_delila_pulser/RAW/SDataR_run_same_settings_delila_pulser.root");

	// TFile *f_positron = TFile::Open("/home/georgen24/Documents/LICENTA/coduri/Root_files/run_same_settings_delila_60Co/RAW/SDataR_run_same_settings_delila_60Co.root");
	TTree *t_data = (TTree *)f_positron->Get("Data_R");
	UShort_t t_channel;
	UInt_t entries;

	entries = t_data->GetEntries();
	TArrayS *samples = nullptr;
	t_data->SetBranchAddress("Samples", &samples);

	t_data->SetBranchAddress("Channel", &t_channel);
	Double_t cfd_delay0[992], cfd_att0[992], cfd0[992], cfd_delay1[992], cfd_att1[992], cfd1[992];

	std::vector<double> sampleNumber;
	std::vector<double> rawSample;
	TGraph *ch0;
	int baseline = 0;
	int iter = 0;
	for (int iEntry = 0; iEntry < 1; iEntry++)
	{

		t_data->GetEntry(iEntry);
		if (t_channel == 1)
		{

			if (samples != nullptr)
			{
				for (int j = 0; j < samples->GetSize() - 1; j++)
				{
					sampleNumber.push_back(j);
					rawSample.push_back(static_cast<double>((samples->At(j) + samples->At(j + 1)) / 2));
				}
			}

			for (int k = 0; k < 100; k++)
			{
				baseline += rawSample[k];
			}
			baseline /= 100;
		}

		std::vector<double> CFD_semnal1 = CFD(rawSample, sampleNumber, baseline, 0.75, 8);
		ch0 = new TGraph(sampleNumber.size(), sampleNumber.data(), CFD_semnal1.data());

		ch0->SetTitle("CFD Signal");
		ch0->Draw("APL*");
		ch0->GetXaxis()->SetTitle("Samples");
		ch0->GetYaxis()->SetTitle("ADC Channels");
		ch0->GetXaxis()->SetRangeUser(200, 260);

		TLine *line = new TLine(200, 0, 260, 0);
		line->SetLineColor(kRed);
		line->SetLineWidth(2);

		line->Draw();
	}
}