/*
 * makeAntSphericalCosineProductsFiles.cxx

 *  Created on: Oct 03, 2020
 *      Author: John Russell
 *
 *  Version of "makeAntSphericalCosineProductsFiles" to be compiled.
 */

using namespace std;

#include <cmath>
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"
//#include "DeltaT.h"
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


void makeAntSphericalCosineProduct(TFile * outFile, AnitaPol::AnitaPol_t pol, int ant1, int ant2) {

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

	TH2D antSphCosProduct(TString::Format("%s_%d_%d", pol == AnitaPol::kHorizontal ? "H" : "V", ant1, ant2), TString::Format("%s_%d_%d", pol == AnitaPol::kHorizontal ? "H" : "V", ant1, ant2), NPhi, minPhi, maxPhi, NNegTheta, minNegTheta, maxNegTheta);

	const UCorrelator::AntennaPositions * ap = UCorrelator::AntennaPositions::instance();

	double phi1 = ap -> phiAnt[pol][ant1] * DEG2RAD;

	double cosTheta1 = cos(M_PI / 18);
	double sinTheta1 = sin(M_PI / 18);

	double phi2 = ap -> phiAnt[pol][ant2] * DEG2RAD;

	double cosTheta2 = cosTheta1;
	double sinTheta2 = sinTheta1;

	//  Fill the TH2D object.
	for (int phiIdx = 1; phiIdx <= antSphCosProduct.GetNbinsX(); ++phiIdx) {

		double binPhi = antSphCosProduct.GetXaxis() -> GetBinCenter(phiIdx) * DEG2RAD;

		double cosDPhi1 = cos(binPhi - phi1);  //  For faster calculation analogous to Abby's method.
		double cosDPhi2 = cos(binPhi - phi2);

		for (int negThetaIdx = 1; negThetaIdx <= antSphCosProduct.GetNbinsY(); ++negThetaIdx) {

			double binNegTheta = antSphCosProduct.GetYaxis() -> GetBinCenter(negThetaIdx) * DEG2RAD;

			double cosTheta = cos(binNegTheta);
			double sinTheta = -sin(binNegTheta);

			double sphCos1 = cosTheta * cosTheta1 * cosDPhi1 + sinTheta * sinTheta1;
			sphCos1 = sphCos1 > 0 ? sphCos1 : 0;  //  Only include positive-valued spherical cosines.
//			sphCos1 /= pol == AnitaPol::kHorizontal ? cosDPhi1 : sinTheta * sinTheta1 * cosDPhi1 + cosTheta * cosTheta1;
//			sphCos1 = sphCos1 > 0 ? sphCos1 : 0;

//			double sphCos1 = pol == AnitaPol::kHorizontal ? cosDPhi1 : sinTheta * sinTheta1 * cosDPhi1 + cosTheta * cosTheta1;
//			sphCos1 = sphCos1 > 0 ? sphCos1 : 0;
//			sphCos1 /= cosTheta * cosTheta1 * cosDPhi1 + sinTheta * sinTheta1;
//			sphCos1 = sphCos1 > 0 ? sphCos1 : 0;

			double sphCos2 = cosTheta * cosTheta2 * cosDPhi2 + sinTheta * sinTheta2;
			sphCos2 = sphCos2 > 0 ? sphCos2 : 0;
//			sphCos2 /= pol == AnitaPol::kHorizontal ? cosDPhi2 : sinTheta * sinTheta2 * cosDPhi2 + cosTheta * cosTheta2;
//			sphCos2 = sphCos2 > 0 ? sphCos2 : 0;

//			double sphCos2 = pol == AnitaPol::kHorizontal ? cosDPhi2 : sinTheta * sinTheta2 * cosDPhi2 + cosTheta * cosTheta2;
//			sphCos2 = sphCos2 > 0 ? sphCos2 : 0;
//			sphCos2 /= cosTheta * cosTheta2 * cosDPhi2 + sinTheta * sinTheta2;
//			sphCos2 = sphCos2 > 0 ? sphCos2 : 0;

			int coverageBinIdx = antSphCosProduct.GetBin(phiIdx, negThetaIdx);
			double coverageBinContent = antSphCosProduct.GetBinContent(coverageBinIdx);

			antSphCosProduct.SetBinContent(coverageBinIdx, coverageBinContent);
			antSphCosProduct.AddBinContent(coverageBinIdx, sphCos1 * sphCos2);
		}
	}

	//  Label the axes.
	antSphCosProduct.GetXaxis() -> SetTitle("#phi");
	antSphCosProduct.GetYaxis() -> SetTitle("-#theta");

	//  Write filled TH2D object to file.
	outFile -> cd();
	antSphCosProduct.Write();
}


int main(int argc, char * argv[]) {

	//  Create file in which coarsely binned TH2D objects will be written.
	TFile coarseFile("antSphericalCosineProductsCoarse.root", "recreate");

	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		for (int j = i; j < NUM_SEAVEYS; ++j) {

			makeAntSphericalCosineProduct(& coarseFile, AnitaPol::kHorizontal, i, j);
			makeAntSphericalCosineProduct(& coarseFile, AnitaPol::kVertical, i, j);
		}

	}

	//  Close the file.
	coarseFile.Close();

	//  Create file in which finely binned TH2D objects will be written.
	TFile fineFile("antSphericalCosineProductsFine.root", "recreate");

	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		for (int j = i; j < NUM_SEAVEYS; ++j) {

			makeAntSphericalCosineProduct(& fineFile, AnitaPol::kHorizontal, i, j);
			makeAntSphericalCosineProduct(& fineFile, AnitaPol::kVertical, i, j);
		}
	}

	//  Close the file.
	fineFile.Close();

	return 1;
}
