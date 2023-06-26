
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1D.h"
#include <TSystem.h>
#include <TROOT.h>
#include "TSpline.h"

std::vector<double> Trapezoid(std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
{
	std::vector<double> result;

	result.push_back(signal_NB[0]);

	for (int i = 1; i < signal_NB.size() - 1; i++)
	{
		double term2 = signal_NB[i - TR_rise_time_K];
		double term3 = signal_NB[i - Trap_Delay_L];
		double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

		result.push_back(result.back() + signal_NB[i] - term2 - term3 + term4);
	}

	return result;
}
// std::vector<double> Trapezoid(std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
// {
// 	std::vector<double> result, p, s, y, d;

// 	result.push_back(signal_NB[0]);
// 	p.push_back(signal_NB[0]);
// 	s.push_back(signal_NB[0]);
// 	y.push_back(signal_NB[0]);

// 	for (int i = 1; i < signal_NB.size(); i++)
// 	{
// 		double term2 = signal_NB[i - TR_rise_time_K];
// 		double term3 = signal_NB[i - Trap_Delay_L];
// 		double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

// 		double di = signal_NB[i] - term2 - term3 + term4 - signal_NB[i - 1] + term2 + term3 - term4;
// 		d.push_back(di);

// 		result.push_back(result.back() + di);
// 		p.push_back(p.back() + di);
// 		s.push_back(s.back() + p[i] + 500 * di);
// 		y.push_back(y.back() + result[i]);
// 	}

// 	return y;
// }
std::vector<double>
Baseline(const std::vector<double> signal)
{
	std::vector<double> signal_NB;
	int baseline = std::accumulate(signal.begin(), signal.begin() + 200, 0.0) / 200;

	std::cout << "Baseline is " << baseline << std::endl;

	// for (int j = 0; j < 2500; j++)
	// {
	// signal_NB.push_back(0);
	// }

	std::vector<double> x;
	for (int i = 0; i < signal.size(); i++)
	{
		signal_NB.push_back(signal[i] - baseline);
	}

	return signal_NB;
}

std::vector<double> Differentiate_Signal(const std::vector<double> signal, double kI, double kD)
{
	std::vector<double> signal_diff, yN;
	signal_diff.clear();
	yN.clear();
	signal_diff.push_back(signal[0]);
	yN.push_back(signal[0]);
	std::vector<double> x;
	// Double_t R = 1000;
	// Double_t C = 1e-6;
	/* double kI = 0.0;
	double kD = 0.0;

	for (auto i = 0.0; i < 0.1; i = i + 0.001)
	{
		kI = i;
		kD = i / 10;
	}; */

	// double kI = 0.007;
	// double kD = 0.0900;

	for (int i = 1; i < signal.size() - 1; i++)
	{

		yN.push_back((yN.back() + signal[i] - signal[i - 1]) / (1 + kD));
		signal_diff.push_back((signal_diff.back() + (1 + kI) * yN[i] - yN[i - 1]) / (1 + kD));

		// signal_diff.push_back((signal_diff.back() + (1 + kI) * signal[i] - signal[i - 1]) / (1 + kD));
		// signal_diff.push_back(signal_diff.back() + (1 + kI) * signal[i] - signal[i - 1]);
	}

	return signal_diff;
}

void RCCR2()
{

	TString filename = "test_000036.root";
	TString treename = "events";
	auto run_file = new TFile(filename, "READ");
	TTree *run_tree;
	TCanvas *c1 = new TCanvas("c1", "Exclusion graphs examples", 200, 10, 600, 400);
	c1->SetGrid();

	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Signals");

	// loc_pad->Divide(3, 1);
	// loc_pad->Draw();

	run_tree = run_file->Get<TTree>(treename);

	Int_t channel;
	ULong_t timestamp;
	Short_t energy;

	vector<double> *waveform1 = {};
	vector<double> *waveform2 = {};
	UInt_t entries, i;

	Double_t TR_Flat_top_M = 500, TR_rise_time_K = 2500;
	double_t Trap_Delay_L = TR_Flat_top_M + TR_rise_time_K;

	run_tree->SetBranchAddress("channel", &channel);
	run_tree->SetBranchAddress("energy", &energy);
	run_tree->SetBranchAddress("timestamp", &timestamp);
	run_tree->SetBranchAddress("waveform1", &waveform1);
	run_tree->SetBranchAddress("waveform2", &waveform2);

	entries = run_tree->GetEntries();
	printf("entries = %d\n", entries);

	std::vector<double> x, raw_signal, trapezoid_ifin, signal_diff;

	for (int i = 5; i < 6; i++)
	{
		run_tree->GetEntry(i);
		cout << "ce canal e " << channel << endl;
		int bin = 0;
		for (int j = 0; j < waveform2->size() - 8; j++)
		{
			x.push_back(j);
			raw_signal.push_back(((*waveform2)[j] + (*waveform2)[j + 1] + (*waveform2)[j + 2] + (*waveform2)[j + 3]) / 4);
			// raw_signal.push_back((*waveform2)[j]);
			// trapezoid_ifin.push_back((*waveform1)[j] / 2);

			trapezoid_ifin.push_back(((*waveform1)[j] + (*waveform1)[j + 1] + (*waveform1)[j + 2] + (*waveform1)[j + 3] + (*waveform1)[j + 4] + (*waveform1)[j + 5] + (*waveform1)[j + 6] + (*waveform1)[j + 7]) / 8);
		}
	}
	std::cout << "test" << std::endl;
	std::vector<double> raw_signal_NB = Baseline(raw_signal);
	std::cout << "test2" << std::endl;
	//	std::vector<double> trapezoid_NB = Baseline(trapezoid_ifin);

	double kI = 0.01;
	double kD = 0.000;

	int k = 0;

	std::vector<double> output_signal_adder;

	std::vector<double> differential_1 = Differentiate_Signal(raw_signal_NB, 0.02, 0.0090);
	std::vector<double> differential_2 = Differentiate_Signal(raw_signal_NB, 0.0121, 0.121);
	std::vector<double> differential_3 = Differentiate_Signal(raw_signal_NB, 0.02, 0.0082);

	std::vector<double> result = Trapezoid(raw_signal_NB, TR_rise_time_K, Trap_Delay_L);

	auto sig_graph_1 = new TGraph(x.size(), x.data(), result.data());

	sig_graph_1->SetMarkerColor(kRed);
	sig_graph_1->SetMarkerStyle(kFullSquare);

	//   sig_graph_3->GetYaxis()->SetRange(-20, 20);
	auto sig_graph_2 = new TGraph(x.size(), x.data(), differential_2.data());
	sig_graph_2->SetMarkerColor(kGreen);
	auto sig_graph_3 = new TGraph(x.size(), x.data(), differential_3.data());
	sig_graph_3->SetMarkerColor(kBlue);
	auto sig_graph_4 = new TGraph(x.size(), x.data(), raw_signal_NB.data());
	sig_graph_4->SetMarkerColor(kMagenta);

	mg->Add(sig_graph_1);
	// mg->Add(sig_graph_2);
	//	mg->Add(sig_graph_3);
	// mg->Add(sig_graph_4);
	mg->Draw("AP*");

	output_signal_adder.clear();
}
