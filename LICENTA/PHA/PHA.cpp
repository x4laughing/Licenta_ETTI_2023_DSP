
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
// #include "header.h"
// std::vector<double> Trapezoid(std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
// {
// std::vector<double> result, p, s, y, d;

// result.push_back(signal_NB[0]);
// p.push_back(signal_NB[0]);
// s.push_back(signal_NB[0]);
// y.push_back(signal_NB[0]);

// for (int i = 1; i < signal_NB.size(); i++)
// {
// double term2 = signal_NB[i - TR_rise_time_K];
// double term3 = signal_NB[i - Trap_Delay_L];
// double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

// double di = signal_NB[i] - term2 - term3 + term4 - signal_NB[i - 1] + term2 + term3 - term4;
// d.push_back(di);

// result.push_back(result.back() + di);
// p.push_back(p.back() + di);
// s.push_back(s.back() + p[i] + 500 * di);
// y.push_back(y.back() + result[i]);
// }

// return y;
// }

// std::vector<double> Trapezoid(std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
void Trapezoid(std::vector<double> *result, std::vector<double> signal_NB, int TR_rise_time_K, int Trap_Delay_L)
{
	// std::vector<double> result;
	result->clear();

	result->push_back(signal_NB[0]);

	for (int i = 0; i < signal_NB.size(); i++)
	{
		double term2 = signal_NB[i - TR_rise_time_K];
		double term3 = signal_NB[i - Trap_Delay_L];
		double term4 = signal_NB[i - Trap_Delay_L - TR_rise_time_K];

		result->push_back(result->back() + signal_NB[i] - term2 - term3 + term4);
	}

	// return result;
}

// std::vector<double> trapezoidal_energy(const std::vector<double> &x, const std::vector<double> &y)
// {
// std::vector<double> energy_spectrum;
// double energy = 0;
// for (int i = 1700; i < 2000; i++)
// {
// double segment_width = x[i + 1] - x[i];
// double segment_height_avg = (y[i + 1] + y[i]) / 2;
// double segment_area = segment_width * segment_height_avg;
// energy += segment_area;
// energy_spectrum.push_back(energy);
// }
// return energy_spectrum;
// }
std::vector<double> Baseline(const std::vector<double> signal)
{
	std::vector<double> signal_NB;
	signal_NB.clear();
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

std::vector<double> Differentiate_Signal(const std::vector<double> signal)
{ // const double kI = 0.007, const double kD = 0.07
	std::vector<double> signal_diff, yN;
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
	// double kD = 0.0090;
	double kI = 0.0121;
	double kD = 0.121;
	for (int i = 1; i < signal.size() - 1; i++)
	{

		yN.push_back((yN.back() + signal[i] - signal[i - 1]) / (1 + kD));
		signal_diff.push_back((signal_diff.back() + (1 + kI) * yN[i] - yN[i - 1]) / (1 + kD));

		// signal_diff.push_back((signal_diff.back() + (1 + kI) * signal[i] - signal[i - 1]) / (1 + kD));
		// signal_diff.push_back(signal_diff.back() + (1 + kI) * signal[i] - signal[i - 1]);
	}

	return signal_diff;
}
TPad *createThreeGraphsOnPad(const std::vector<double> &x, const std::vector<double> &first_signal, const std::vector<double> &second_signal, const std::vector<double> &third_signal)
{
	auto pad = new TPad("pad name", "pad title", 0, 0, 1, 1);
	pad->Divide(3, 1, 0.01, 0.01);
	pad->Draw();
	pad->cd(1);

	auto sig_graph_1 = new TGraph(x.size(), x.data(), first_signal.data());
	pad->cd(1);
	sig_graph_1->SetTitle(Form("Differential signal"));
	sig_graph_1->Draw("APL");

	auto sig_graph_2 = new TGraph(x.size(), x.data(), second_signal.data());
	sig_graph_2->SetMarkerColor(kBlue);
	sig_graph_2->SetMarkerStyle(kFullCircle);
	sig_graph_2->SetTitle(Form("Tapezoid PZC "));
	// sig_graph_2->GetXaxis()->SetRange(2400, 3100);
	// sig_graph_2->GetYaxis()->SetRange(2595, 2615);
	pad->cd(2);
	sig_graph_2->Draw("AP*");

	auto sig_graph_3 = new TGraph(x.size(), x.data(), third_signal.data());
	sig_graph_3->SetMarkerColor(kRed);
	sig_graph_3->SetMarkerStyle(kFullSquare);
	sig_graph_3->SetTitle(Form("Tapezoid RAW"));
	pad->cd(3);
	sig_graph_3->Draw("AP*");

	return pad;
}

void PHA()
{

	TString filename = "test_000031.root";
	TString treename = "events";
	auto run_file = new TFile(filename, "READ");
	TTree *run_tree;
	auto loc_canv = new TCanvas("runs_from_DAQ.pdf", "c1");
	auto loc_pad = new TPad("pad name", "pad title", 0, 0, 1, 1);
	loc_pad->Divide(3, 1);
	loc_pad->Draw();

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

	std::vector<double> x, raw_signal, trapezoid_ifin, analog_adder, signal_diff2, output_signal_adder_negative, output_signal_adder, result;

	for (int i = 1; i < 2; i++)
	{
		run_tree->GetEntry(i);
		cout << "ce canal e " << channel << endl;
		int bin = 0;
		for (int j = 0; j < waveform2->size() - 8; j++)
		{
			x.push_back(j);
			raw_signal.push_back(((*waveform2)[j] + (*waveform2)[j + 1] + (*waveform2)[j + 2] + (*waveform2)[j + 3]) / 4);
			// raw_signal.push_back((*waveform2)[j]);
			// trapezoid_ifin.push_back((*waveform1)[j]);

			trapezoid_ifin.push_back(((*waveform1)[j] + (*waveform1)[j + 1] + (*waveform1)[j + 2] + (*waveform1)[j + 3] + (*waveform1)[j + 4] + (*waveform1)[j + 5] + (*waveform1)[j + 6] + (*waveform1)[j + 7]) / 8);
		}
	}

	std::vector<double> raw_signal_NB = Baseline(raw_signal);

	// std::vector<double> result = Trapezoid(raw_signal_NB, TR_rise_time_K, Trap_Delay_L);

	std::vector<double> signal_diff = Differentiate_Signal(raw_signal_NB);

	for (int i = 0; i < raw_signal_NB.size(); ++i)
	{
		signal_diff2.push_back(raw_signal_NB[i] * (0.019));
		// signal_diff2.push_back(raw_signal_NB[i] * (1));
	}

	for (int i = 0; i < signal_diff.size(); ++i)
	{
		analog_adder.push_back(signal_diff2[i] + signal_diff[i]);
	}
	// std::vector<double> result = Trapezoid(signal_diff, TR_rise_time_K, Trap_Delay_L);

	// Double_t baselinesignal_diff;
	// for (int i = 2000; i < signal_diff.size(); ++i)
	// {
	// 	baselinesignal_diff = (signal_diff[i] + signal_diff[i + 1]) / 2;
	// 	// cout << "Valoare semnal signal diff e" << baselinesignal_diff << std::endl;
	// }
	// cout << "Valoare semnal signal diff e " << baselinesignal_diff << std::endl;

	// Double_t baselinanalog_adder;
	// for (int i = 2000; i < analog_adder.size(); ++i)
	// {
	// 	baselinanalog_adder = (analog_adder[i] + analog_adder[i + 1]) / 2;
	// 	cout << "Valoare semnal signal diff e" << baselinesignal_diff << std::endl;
	// }
	// cout << "Valoare semnal baselinanalog_addere " << baselinanalog_adder << std::endl;
	// std::vector<double> output_signal_adder;
	Trapezoid(&output_signal_adder, analog_adder, TR_rise_time_K, Trap_Delay_L);
	Trapezoid(&result, raw_signal_NB, TR_rise_time_K, Trap_Delay_L);
	// std::vector<double> result_2 = Trapezoid(output_signal_adder, TR_rise_time_K, Trap_Delay_L);
	for (int i = 0; i < signal_diff.size(); ++i)
	{
		output_signal_adder_negative.push_back(output_signal_adder[i] * (-1.0));
	}
	std::cout << "check1 " << std::endl;

	std::vector<double> final_result, trap_ifin_Nb;

	static constexpr double amplitudine = 7142.85714;
	for (int i = 0; i < output_signal_adder_negative.size(); ++i)
	{
		final_result.push_back(output_signal_adder_negative[i] / (amplitudine));
	}

	static constexpr double BASELINE_IFIN = 216.0; // DE CALCULAT PT FIEACRE
	for (int i = 0; i < trapezoid_ifin.size(); ++i)
	{
		trap_ifin_Nb.push_back(trapezoid_ifin[i] - (BASELINE_IFIN));
	}

	// loc_pad = createThreeGraphsOnPad(x, signal_diff, trapezoid_ifin, output_signal_adder);

	// std::vector<double> signal_diff = Differentiate_Signal(raw_signal_NB);
	// std::vector<double> signal_diff = Baseline(signal_diff);

	// std::vector<double> trapezoidal_spectrum = trapezoidal_energy(x, trapezoid_ifin);
	// std::vector<double> trapezoidal_spectrum_1 = trapezoidal_energy(x, output_signal_adder_negative);

	TPad *pad = createThreeGraphsOnPad(x, signal_diff, trap_ifin_Nb, result);
};

// void PHA_git()
// {
// 	int iter = 0;
// 	while (1)
// 	{
// 		PHA2_git();
// 		iter = iter + 1;
// 		if (iter > 1)
// 			break;
// 	}
// }
// flat top 1 us-> 1000ns->500/4 samples
// traprise 5 us-> 5000 ns -> 2500/4 samples
// peaking 80%-> 800ns -> 400/4 samples
// holdoff 496
// peak holdoff 0.960 us
// rec 20000 ns

// TOP-> m = 500 samples
// BOTTOM -> 2k+m=

// SE IMPARTE la 2