/*
 * makeDeltaTMapFiles.cxx

 *  Created on: Oct 03, 2020
 *      Author: John Russell
 *
 *  Version of "makeDeltaTMapFiles" to be compiled.
 */

using namespace std;

#include <cmath>
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"
#include "DeltaT.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"


#ifndef DEG2RAD
#define DEG2RAD M_PI / 180
#endif

#ifndef RAD2DEG
#define RAD2DEG 180 / M_PI
#endif

const int NPhiCoarse = 180;
const double minPhiCoarse = 0;
const double maxPhiCoarse = 360;
const int NNegThetaCoarse = 100;
const double minNegThetaCoarse = -60;
const double maxNegThetaCoarse = 40;

const int NPhiFine = 760;
const int minPhiFine = -10;
const int maxPhiFine = 370;
const int NNegThetaFine = 240;
const double minNegThetaFine = -70;
const double maxNegThetaFine = 50;


void makeDeltaTMap(TFile * outFile, AnitaPol::AnitaPol_t pol, int ant1, int ant2) {

	//  Don't correlate antennas more than 2 phi sectors apart.
	int phiSep = abs(ant1 - ant2) % NUM_PHI;
	phiSep = min(phiSep, NUM_PHI - phiSep);
	if (phiSep > 2) return;

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	TString outFileNameStr = TString(outFile -> GetName());
	bool isFine = outFileNameStr.Contains("Fine", TString::kIgnoreCase);

	int NPhi = isFine ? NPhiFine : NPhiCoarse;
	double minPhi = isFine ? minPhiFine : minPhiCoarse;
	double maxPhi = isFine ? maxPhiFine : maxPhiCoarse;
	int NNegTheta = isFine ? NNegThetaFine : NNegThetaCoarse;
	double minNegTheta = isFine ? minNegThetaFine : minNegThetaCoarse;
	double maxNegTheta = isFine ? maxNegThetaFine : maxNegThetaCoarse;

	TH2D deltaTMap(TString::Format("%s_%d_%d", pol == AnitaPol::kHorizontal ? "H" : "V", ant1, ant2), TString::Format("%s_%d_%d", pol == AnitaPol::kHorizontal ? "H" : "V", ant1, ant2), NPhi, minPhi, maxPhi, NNegTheta, minNegTheta, maxNegTheta);

	//  Fill the TH2D object.
	for (int phiIdx = 1; phiIdx <= deltaTMap.GetNbinsX(); ++phiIdx) {

		double binPhi = deltaTMap.GetXaxis() -> GetBinCenter(phiIdx) * DEG2RAD;

		for (int negThetaIdx = 1; negThetaIdx <= deltaTMap.GetNbinsY(); ++negThetaIdx) {

			double binNegTheta = deltaTMap.GetYaxis() -> GetBinCenter(negThetaIdx) * DEG2RAD;

			double deltaT = UCorrelator::getDeltaT(ant1, ant2, binPhi * RAD2DEG, -binNegTheta * RAD2DEG, pol);  // Geometric time delay.

			int coverageBinIdx = deltaTMap.GetBin(phiIdx, negThetaIdx);
			double coverageBinContent = deltaTMap.GetBinContent(coverageBinIdx);

			deltaTMap.SetBinContent(coverageBinIdx, coverageBinContent);
			deltaTMap.AddBinContent(coverageBinIdx, deltaT);
		}
	}

	//  Label the axes.
	deltaTMap.GetXaxis() -> SetTitle("#phi");
	deltaTMap.GetYaxis() -> SetTitle("-#theta");

	//  Write filled TH2D object to file.
	outFile -> cd();
	deltaTMap.Write();
}


int main(int argc, char * argv[]) {

	//  Create file in which coarsely binned TH2D objects will be written.
	TFile coarseFile("deltaTMapCoarse.root", "recreate");

	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		for (int j = i + 1; j < NUM_SEAVEYS; ++j) {

			makeDeltaTMap(& coarseFile, AnitaPol::kHorizontal, i, j);
			makeDeltaTMap(& coarseFile, AnitaPol::kVertical, i, j);
		}

	}

	//  Close the file.
	coarseFile.Close();

	//  Create file in which finely binned TH2D objects will be written.
	TFile fineFile("deltaTMapFine.root", "recreate");

	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		for (int j = i + 1; j < NUM_SEAVEYS; ++j) {

			makeDeltaTMap(& fineFile, AnitaPol::kHorizontal, i, j);
			makeDeltaTMap(& fineFile, AnitaPol::kVertical, i, j);
		}
	}

	//  Close the file.
	fineFile.Close();

	return 1;
}
