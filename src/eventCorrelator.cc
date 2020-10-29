/*
 * eventCorrelator.cc
 *
 *  Created on: Jul 21, 2020
 *      Author: John Russell
 *
 *  Header file from which to reference correlator functions.
 */

#include "eventCorrelator.h"

#ifndef DEG2RAD
#define DEG2RAD M_PI / 180
#endif

#ifndef RAD2DEG
#define RAD2DEG 180 / M_PI
#endif

#ifndef PHI_SECTOR_ANGLE
#define PHI_SECTOR_ANGLE 360 / NUM_PHI
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

const int NPhiZoom = 40;
const double dPhiZoom = 0.5;
const int NNegThetaZoom = 40;
const double dNegThetaZoom = 0.5;


TFile deltaTMapCoarseFile("deltaTMapCoarse.root");  //  Because the files are opened outside of functions below, they aren't ever closed. As such, pointers to objects within them must be explicitly deleted when no longer in use, otherwise they will remain in memory.
TFile deltaTMapFineFile("deltaTMapFine.root");

TFile antSphericalCosineProductsCoarseFile("antSphericalCosineProductsCoarse.root");
TFile antSphericalCosineProductsFineFile("antSphericalCosineProductsFine.root");


/*
 * Get antenna pairs used to construct interferometric maps.
 */
vector<vector<int>> getAntennaPairs(vector<int> neighboringAntennas) {

	//  Check if all neighboring antennas allowed.
	if (neighboringAntennas.empty()) {

		neighboringAntennas.resize(NUM_SEAVEYS);
		for (int i = 0; i < NUM_SEAVEYS; ++i) neighboringAntennas[i] = i;
	}

	vector<vector<int>> antennaPairs;  //  Where antenna pairs will be placed.

	//  Nested for loop which fills the antennaPairs vector.
	for (int i = 0; i < neighboringAntennas.size(); ++i) {

		for (int j = i + 1; j < neighboringAntennas.size(); ++j) {

			int phiSep = abs(neighboringAntennas[i] - neighboringAntennas[j]) % NUM_PHI;
			phiSep = min(phiSep, NUM_PHI - phiSep);
			if (phiSep > 2) continue;  //  Exclude pairs more than 2 phi sectors apart.

			antennaPairs.push_back({neighboringAntennas[i], neighboringAntennas[j]});
		}
	}

	return antennaPairs;
}


/*
 * Given an event, cons	antSphericalCosineProductsFile.Close();
 * truct the cross-correlation between two antennas for a given polarization.
 */
void fillMapsPair(FilteredAnitaEvent * filtEvent, TH2D * responseMap, int ant1, int ant2, AnitaPol::AnitaPol_t pol) {

	//  Check dimensions of responseMap, to determine which spherical cosine pair file to reference.
	int NPhi = responseMap -> GetNbinsX();
	int NNegTheta = responseMap -> GetNbinsY();

	//  Char associated with polarization.
	const char * polChar = (pol == AnitaPol::kHorizontal) ? "H" : "V";

	//  Determine which pair of spherical cosines histogram to reference.
		TH2D * deltaTMap = 0;
		if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

			deltaTMap = (TH2D *) deltaTMapCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));
	//		deltaTMap -> Scale(-1);  //  Testing reciprocity. Normally this is commented out.

			if (!deltaTMap) {  //  Check if NULL pointer.

				deltaTMap = (TH2D *) deltaTMapCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
				deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
			}

		} else {

			deltaTMap = (TH2D *) deltaTMapFineFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));

			if (!deltaTMap) {

				deltaTMap = (TH2D *) deltaTMapFineFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
				deltaTMap -> Scale(-1);
			}
		}

		TH2D * antSphCosProduct = 0;
		if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

			antSphCosProduct = (TH2D *) antSphericalCosineProductsCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));

			if (!antSphCosProduct) antSphCosProduct = (TH2D *) antSphericalCosineProductsCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));

		} else {

			antSphCosProduct = (TH2D *) antSphericalCosineProductsFineFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));

			if (!antSphCosProduct) antSphCosProduct = (TH2D *) antSphericalCosineProductsFineFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
		}

//	//  Accessing relevant TH2D objects.
//	TFile deltaTMapFile = (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? deltaTMapCoarseFile : deltaTMapFineFile;
//	TH2D * deltaTMap = (TH2D *) deltaTMapFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));
//	if (!deltaTMap) {
//
//		deltaTMap = (TH2D *) deltaTMapFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
//		deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
//	}
////	deltaTMap -> SetDirectory(0);
//
//	TFile antSphericalCosineProductsFile= (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? antSphericalCosineProductsCoarseFile : antSphericalCosineProductsFineFile;
//	TH2D * antSphCosProduct = (TH2D *) antSphericalCosineProductsFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));
//	if (!antSphCosProduct) antSphCosProduct = (TH2D *) antSphericalCosineProductsFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
////	antSphCosProduct -> SetDirectory(0);

	//  Get index differences between reseponseMap and antSphCosProduct. Relevant for zoomed maps.
	int dPhiIdx = int((responseMap -> GetXaxis() -> GetBinCenter(1) - deltaTMap -> GetXaxis() -> GetBinCenter(1)) / dPhiZoom);
	int dNegThetaIdx = int((responseMap -> GetYaxis() -> GetBinCenter(1) - deltaTMap -> GetYaxis() -> GetBinCenter(1)) / dNegThetaZoom);

	//  Setting up waveforms and correlations to be used in interferometric map.
	TGraph waveformAnt1 = TGraph(* filtEvent -> getFilteredGraph(ant1, pol) -> even());
	TGraph waveformAnt2 = TGraph(* filtEvent -> getFilteredGraph(ant2, pol) -> even());

	TGraph crossCorr = getCorrGraph(& waveformAnt1, & waveformAnt2);
	crossCorr.SetBit(TGraph::kIsSortedX);  //  This should significantly expedite interpolation.
	TSpline3 crossCorrSpline("crossCorrSpline", & crossCorr);  //  Expedite processing by evaluating the spline once, then referencing it.

	for (int phiIdx = 1; phiIdx <= NPhi; ++phiIdx) {

		double binPhi = responseMap -> GetXaxis() -> GetBinCenter(phiIdx) * DEG2RAD;

		for (int negThetaIdx = 1; negThetaIdx <= NNegTheta; ++negThetaIdx) {

			double rho = antSphCosProduct -> GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx);

			if (rho <= 0) continue;  //  Reference the spherical cosine pair to determine if bin should be filled.

			double binNegTheta = responseMap -> GetYaxis() -> GetBinCenter(negThetaIdx) * DEG2RAD;

			double deltaT = deltaTMap ? deltaTMap -> GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx) : 0;

 			double response = crossCorrSpline.Eval(deltaT);

			int responseBinIdx = responseMap -> GetBin(phiIdx, negThetaIdx);
			double responseBinContent = responseMap -> GetBinContent(responseBinIdx);
			responseMap -> SetBinContent(responseBinIdx, responseBinContent);

			if (ant1 != ant2) responseMap -> AddBinContent(responseBinIdx, rho * response);
			else responseMap -> AddBinContent(responseBinIdx, rho * response / 2);
		}
	}

//	//  Delete pointers to file objects, then close the files.
//	gROOT -> cd();
	delete deltaTMap;
	delete antSphCosProduct;
//
//	gROOT -> cd();
//	deltaTMapFile.Close();
//	antSphericalCosineProductsFile.Close();
}


/*
 * Function which iteratively calls "fillMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillMapsAll(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<vector<int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	responseMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillMapsPair(filtEvent, & mapPair[i], antennaPairs[i][0], antennaPairs[i][1], pol);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) responseMap -> Add(& mapPair[i]);

	responseMap -> GetXaxis() -> SetTitle("#phi");
	responseMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * For a given event, return a 2-dimensional vector containing sums of square voltages for each antenna. Vector has dimensions of
 * [2][NUM_SEAVEYS]. where first index corresponds to 0 for Hpol, 1 for Vpol.
 */
vector<vector<double>> getEventTotalPowers(FilteredAnitaEvent * filtEvent) {

	vector<vector<double>> totalPower(2, vector<double>(NUM_SEAVEYS));

	#pragma omp parallel for
	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		TGraph waveformAntH = TGraph(* filtEvent -> getFilteredGraph(i, AnitaPol::kHorizontal) -> even());
		TGraph waveformAntV = TGraph(* filtEvent -> getFilteredGraph(i, AnitaPol::kVertical) -> even());

		for (int j = 0; j < waveformAntH.GetN(); ++j) totalPower[0][i] += waveformAntH.GetY()[j] * waveformAntH.GetY()[j];
		for (int j = 0; j < waveformAntV.GetN(); ++j) totalPower[1][i] += waveformAntV.GetY()[j] * waveformAntV.GetY()[j];
	}

	return totalPower;
}


/*
 * For normalization purposes. For both H and Vpol channel of each antenna, calculate the "total power" histogram.
 * For a pair of antennas at a given azimuth and elevation, normalization is given by taking the product of the two "total power" histograms
 * at that point and taking the square root of that product.
 */
vector<TH2D> getTotalPowerMaps(FilteredAnitaEvent * filtEvent, bool isFine) {

	//  Get total powers to be used.
	vector<vector<double>> totalPowers = getEventTotalPowers(filtEvent);

	vector<TH2D> totalPowerMaps(NUM_SEAVEYS);

	//  Get dimensions of coverage map.
	int NPhi = isFine ? NPhiFine : NPhiCoarse;
	double minPhi = isFine ? minPhiFine : minPhiCoarse;
	double maxPhi = isFine ? maxPhiFine : maxPhiCoarse;
	int NNegTheta = isFine ? NNegThetaFine : NNegThetaCoarse;
	double minNegTheta = isFine ? minNegThetaFine : minNegThetaCoarse;
	double maxNegTheta = isFine ? maxNegThetaFine : maxNegThetaCoarse;

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

//	//  Accessing file from which to read histograms.
//	TFile antSphericalCosineProductsFile = isFine ? antSphericalCosineProductsFineFile : antSphericalCosineProductsCoarseFile;

	#pragma omp parallel for
	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		TH2D * hHist = 0, * vHist = 0;

		if (isFine) {

			hHist = (TH2D *) antSphericalCosineProductsFineFile.Get(TString::Format("H_%d_%d", i, i));
			vHist = (TH2D *) antSphericalCosineProductsFineFile.Get(TString::Format("V_%d_%d", i, i));

		} else {

			hHist = (TH2D *) antSphericalCosineProductsCoarseFile.Get(TString::Format("H_%d_%d", i, i));
			vHist = (TH2D *) antSphericalCosineProductsCoarseFile.Get(TString::Format("V_%d_%d", i, i));
		}

//		TH2D * hHist = (TH2D *) antSphericalCosineProductsFile.Get(TString::Format("H_%d_%d", i, i));
////		hHist -> SetDirectory(0);
//
//		TH2D * vHist = (TH2D *) antSphericalCosineProductsFile.Get(TString::Format("V_%d_%d", i, i));
////		vHist -> SetDirectory(0);

		totalPowerMaps[i] = TH2D(TString::Format("Ant_%d", i), TString::Format("Ant_%d", i), NPhi, minPhi, maxPhi, NNegTheta, minNegTheta, maxNegTheta);
		totalPowerMaps[i].Add(hHist, vHist, totalPowers[0][i], totalPowers[1][i]);

		totalPowerMaps[i].GetXaxis() -> SetTitle("#phi");
		totalPowerMaps[i].GetYaxis() -> SetTitle("-#theta");

//		gROOT -> cd();
		delete hHist;
		delete vHist;
	}

//	//  Close the file.
//	gROOT -> cd();
//	antSphericalCosineProductsFile.Close();

	return totalPowerMaps;
}


/*
 * Each antenna cross-correlation has an upper bound given by the Cauchy-Schwarz inequality. The upper bound here is over both horizontal and vertical polarization.
 */
void fillCoverageMapsPair(TH2D * coverageMap, vector<TH2D> totalPowerMaps, int ant1, int ant2) {

	//  Access coverage maps to use.
	TH2D totalPowerMap1 = totalPowerMaps[ant1];
	TH2D totalPowerMap2 = totalPowerMaps[ant2];

	//  Get index differences between coverageMap and totalPowerMaps. Relevant for zoomed maps.
	int dPhiIdx = int((coverageMap -> GetXaxis() -> GetBinCenter(1) - totalPowerMap1.GetXaxis() -> GetBinCenter(1)) / dPhiZoom);
	int dNegThetaIdx = int((coverageMap -> GetYaxis() -> GetBinCenter(1) - totalPowerMap1.GetYaxis() -> GetBinCenter(1)) / dNegThetaZoom);

	for (int phiIdx = 1; phiIdx <= coverageMap -> GetNbinsX(); ++phiIdx) {

		for (int negThetaIdx = 1; negThetaIdx <= coverageMap -> GetNbinsY(); ++negThetaIdx) {

			double powerVal1 = totalPowerMap1.GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx);
			if (powerVal1 <= 0) continue;

			double powerVal2 = totalPowerMap2.GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx);
			if (powerVal2 <= 0) continue;

			double rho = sqrt(powerVal1 * powerVal2);

			int coverageBinIdx = coverageMap -> GetBin(phiIdx, negThetaIdx);
			double coverageBinContent = coverageMap -> GetBinContent(coverageBinIdx);
			coverageMap -> SetBinContent(coverageBinIdx, coverageBinContent);

			if (ant1 != ant2) coverageMap -> AddBinContent(coverageBinIdx, rho);
			else coverageMap -> AddBinContent(coverageBinIdx, rho / 2);
		}
	}
}


/*
 * Function which iteratively calls "fillCoverageMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillCoverageMapsAll(TH2D * coverageMap, vector<TH2D> totalPowerMaps, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<vector<int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	coverageMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair coverage histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillCoverageMapsPair(& mapPair[i], totalPowerMaps, antennaPairs[i][0], antennaPairs[i][1]);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) coverageMap -> Add(& mapPair[i]);

	coverageMap -> GetXaxis() -> SetTitle("#phi");
	coverageMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * For a given event, calls to "fillMapsAll()" to construct corresponding unnormalized interferoemtric maps.
 */
vector<TH2D> makeUnnormalizedEventInterferometricMaps(int eventNum, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Setting up histograms.
	TH2D responseMapHPol("responseMapHPol", "HPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);
	TH2D responseMapVPol("responseMapVPol", "VPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);
	TH2D responseMapUnpol("responseMapUnpol", "Unpolarized response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);

	//  Running functions with produced histograms.
	fillMapsAll(& filtEvent, & responseMapHPol, AnitaPol::kHorizontal);
	fillMapsAll(& filtEvent, & responseMapVPol, AnitaPol::kVertical);

	responseMapUnpol.Add(& responseMapHPol);
	responseMapUnpol.Add(& responseMapVPol);
	responseMapUnpol.Scale(0.5);
	responseMapUnpol.GetXaxis() -> SetTitle("#phi");
	responseMapUnpol.GetYaxis() -> SetTitle("-#theta");

	//  Produce the vector of TH2D's, then return it.
	vector<TH2D> eventInterferometricMaps(3);
	eventInterferometricMaps[0] = responseMapHPol;
	eventInterferometricMaps[1] = responseMapVPol;
	eventInterferometricMaps[2] = responseMapUnpol;

	return eventInterferometricMaps;
}


/*
 * For a given event, calls to "fillMapsAll()" and "fillCoverageMapsAll()" to construct corresponding interferometric maps.
 */
vector<TH2D> makeEventInterferometricMaps(int eventNum, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Setting up histograms.
	TH2D responseMapHPol("responseMapHPol", "HPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);
	TH2D responseMapVPol("responseMapVPol", "VPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);

	TH2D coverageMap("coverageMap", "Unpolarized antenna pair coverage", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);

	TH2D normalizedResponseMapHPol("normalizedResponseMapHPol", "Normalized HPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);
	TH2D normalizedResponseMapVPol("normalizedResponseMapVPol", "Normalized VPol response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);
	TH2D normalizedResponseMapUnpol("normalizedResponseMapUnpol", "Normalized Unpolarized response", NPhiCoarse, minPhiCoarse, maxPhiCoarse, NNegThetaCoarse, minNegThetaCoarse, maxNegThetaCoarse);

	//  Running functions with produced histograms.
	fillMapsAll(& filtEvent, & responseMapHPol, AnitaPol::kHorizontal);
	fillMapsAll(& filtEvent, & responseMapVPol, AnitaPol::kVertical);

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent);

	fillCoverageMapsAll(& coverageMap, totalPowerMaps);

	normalizedResponseMapHPol.Add(& responseMapHPol);
	normalizedResponseMapHPol.Divide(& coverageMap);
	normalizedResponseMapHPol.Scale(2);
	normalizedResponseMapHPol.GetXaxis() -> SetTitle("#phi");
	normalizedResponseMapHPol.GetYaxis() -> SetTitle("-#theta");

	normalizedResponseMapVPol.Add(& responseMapVPol);
	normalizedResponseMapVPol.Divide(& coverageMap);
	normalizedResponseMapVPol.Scale(2);
	normalizedResponseMapVPol.GetXaxis() -> SetTitle("#phi");
	normalizedResponseMapVPol.GetYaxis() -> SetTitle("-#theta");

	normalizedResponseMapUnpol.Add(& normalizedResponseMapHPol);
	normalizedResponseMapUnpol.Add(& normalizedResponseMapVPol);
	normalizedResponseMapUnpol.Scale(0.5);
	normalizedResponseMapUnpol.GetXaxis() -> SetTitle("#phi");
	normalizedResponseMapUnpol.GetYaxis() -> SetTitle("-#theta");

	//  Produce the vector of TH2D's, then return it.
	vector<TH2D> eventInterferometricMaps(3);
	eventInterferometricMaps[0] = normalizedResponseMapHPol;
	eventInterferometricMaps[1] = normalizedResponseMapVPol;
	eventInterferometricMaps[2] = normalizedResponseMapUnpol;

	return eventInterferometricMaps;
}


/*
 * For a given interferometric map, determine the locations in phi and theta of the peak in the map.
 */
vector<double> getInterferometricPeakLocation(TH2D * responseMap) {

	//  An array returning location in phi (0) and negative theta (1) of that map.
	vector<double> peakLocation(2);

	//  First, determine location of the maximum peak in the interferometric map.
	int peakIdx = responseMap -> GetMaximumBin();

	int peakPhiIdx, peakNegThetaIdx, peakZIdx;

	responseMap -> GetBinXYZ(peakIdx, peakPhiIdx, peakNegThetaIdx, peakZIdx);

	peakLocation[0] = responseMap -> GetXaxis() -> GetBinCenter(peakPhiIdx);
	peakLocation[1] = responseMap -> GetYaxis() -> GetBinCenter(peakNegThetaIdx);

	return peakLocation;
}


/*
 * Given where an interferometric peak occurs, determine which phi sector it is, then return antenna numbers up to two phi sectors apart.
 */
vector<int> getNeighboringAntennas(double peakPhi) {

	//  First, determine which phi sector the peak phi corresponds. Center of first phi sector occurs at -45 degrees (315 degrees),
	//  so it starts at -45-11.25 = -2.5 * PHI_SECTOR_ANGLE = -56.25 degrees (303.75 degrees). Starting with index 0 instead of 1.
	int phiSector = int(FFTtools::wrap(peakPhi + 2.5 * PHI_SECTOR_ANGLE, 360) / PHI_SECTOR_ANGLE);

	vector<int> neighboringAntennas(15);

	int antNum = 0;
	for (int i = 0; i < NUM_PHI; ++i) {

		int phiSep = abs(i - phiSector) % NUM_PHI;
		phiSep = min(phiSep, NUM_PHI - phiSep);

		if (phiSep <= 2) {

			neighboringAntennas[antNum] = i;
			neighboringAntennas[antNum + 5] = i + NUM_PHI;
			neighboringAntennas[antNum + 10] = i + 2 * NUM_PHI;

			++antNum;
		}
	}

	return neighboringAntennas;
}


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 * No normalization applied.
 */
TH2D makePeakUnnormalizedInterferometricMap(TH2D * responseMap, int eventNum, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  First, determine location of the maximum peak in the interferometric map.
	vector<double> peakLocation = getInterferometricPeakLocation(responseMap);

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakLocation[0]);

	//  Set dimensions of interferometric maps.
	double minPhi = peakLocation[0] - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakLocation[0] + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakLocation[1] - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakLocation[1] + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	TString responseMapName = TString(responseMap -> GetName());
	const char * polType;
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) polType = "HPol";
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) polType = "VPol";
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakResponseMap("peakResponseMap", "Response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakResponseMap;
}


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakUnnormalizedInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakPhi);

	//  Set dimensions of interferometric maps.
	double minPhi = peakPhi - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakPhi + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakNegTheta - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakNegTheta + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check polarization type of peak interferometric map to be generated.
	const char * polType;
	if (pol == AnitaPol::kHorizontal) polType = "HPol";
	else if (pol == AnitaPol::kVertical) polType = "VPol";
	else if (pol == AnitaPol::kNotAPol) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakResponseMap("peakResponseMap", "Response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillMapsAll(& filtEvent, & peakResponseMap, pol, neighboringAntennas);
	else {

		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakResponseMap;
}

/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 */
TH2D makePeakInterferometricMap(TH2D * responseMap, int eventNum, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  First, determine location of the maximum peak in the interferometric map.
	vector<double> peakLocation = getInterferometricPeakLocation(responseMap);

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakLocation[0]);

	//  Set dimensions of interferometric maps.
	double minPhi = peakLocation[0] - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakLocation[0] + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakLocation[1] - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakLocation[1] + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	TString responseMapName = TString(responseMap -> GetName());
	const char * polType;
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) polType = "HPol";
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) polType = "VPol";
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakResponseMap("peakResponseMap", "Response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakCoverageMap("peakCoverageMap", "Antenna pair coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakNormalizedResponseMap("peakNormalizedResponseMap", TString::Format("%s normalized response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent, true);

	fillCoverageMapsAll(& peakCoverageMap, totalPowerMaps, neighboringAntennas);

	peakNormalizedResponseMap.Add(& peakResponseMap);
	peakNormalizedResponseMap.Divide(& peakCoverageMap);
	if (!responseMapName.Contains("Unpol", TString::kIgnoreCase)) peakNormalizedResponseMap.Scale(2);

	peakNormalizedResponseMap.GetXaxis() -> SetTitle("#phi");
	peakNormalizedResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakNormalizedResponseMap;
}


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, bool useBroadband, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.
	if (useBroadband) fStrat.addOperation(new IFFTDiffFilter);

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakPhi);

	//  Set dimensions of interferometric maps.
	double minPhi = peakPhi - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakPhi + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakNegTheta - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakNegTheta + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	const char * polType;
	if (pol == AnitaPol::kHorizontal) polType = "HPol";
	else if (pol == AnitaPol::kVertical) polType = "VPol";
	else if (pol == AnitaPol::kNotAPol) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakResponseMap("peakResponseMap", "Response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakCoverageMap("peakCoverageMap", "Antenna pair coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakNormalizedResponseMap("peakNormalizedResponseMap", TString::Format("%s normalized response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillMapsAll(& filtEvent, & peakResponseMap, pol, neighboringAntennas);
	else {

		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillMapsAll(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent, true);

	fillCoverageMapsAll(& peakCoverageMap, totalPowerMaps, neighboringAntennas);

	peakNormalizedResponseMap.Add(& peakResponseMap);
	peakNormalizedResponseMap.Divide(& peakCoverageMap);
	if (pol != AnitaPol::kNotAPol) peakNormalizedResponseMap.Scale(2);

	peakNormalizedResponseMap.GetXaxis() -> SetTitle("#phi");
	peakNormalizedResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakNormalizedResponseMap;
}


/*
 * For purposes of testing noncoplanar baselines, create cross-correlation
 * interferometric maps without curvature correction.
 */
void fillFlatMapsPair(FilteredAnitaEvent * filtEvent, TH2D * responseMap, int ant1, int ant2, AnitaPol::AnitaPol_t pol) {

	//  Check dimensions of responseMap, to determine which spherical cosine pair file to reference.
	int NPhi = responseMap -> GetNbinsX();
	int NNegTheta = responseMap -> GetNbinsY();

	//  Char associated with polarization.
	const char * polChar = (pol == AnitaPol::kHorizontal) ? "H" : "V";

	//  Determine which pair of spherical cosines histogram to reference.
	TH2D * deltaTMap = 0;
	if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

		deltaTMap = (TH2D *) deltaTMapCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));

		if (!deltaTMap) {  //  Check if NULL pointer.

			deltaTMap = (TH2D *) deltaTMapCoarseFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
			deltaTMap -> Scale(-1);  //  Ordering is antisymmetric here.
		}

	} else {

		deltaTMap = (TH2D *) deltaTMapFineFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));

		if (!deltaTMap) {

			deltaTMap = (TH2D *) deltaTMapFineFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
			deltaTMap -> Scale(-1);
		}
	}

//	TFile deltaTMapFile = (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? deltaTMapCoarseFile : deltaTMapFineFile;
//	TH2D * deltaTMap = (TH2D *) deltaTMapFile.Get(TString::Format("%s_%d_%d", polChar, ant1, ant2));
//	if (!deltaTMap) {
//
//		deltaTMap = (TH2D *) deltaTMapFile.Get(TString::Format("%s_%d_%d", polChar, ant2, ant1));
//		deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
//	}
//	deltaTMap -> SetDirectory(0);

	//  Get index differences between reseponseMap and antSphCosProduct. Relevant for zoomed maps.
	int dPhiIdx = int((responseMap -> GetXaxis() -> GetBinCenter(1) - deltaTMap -> GetXaxis() -> GetBinCenter(1)) / dPhiZoom);
	int dNegThetaIdx = int((responseMap -> GetYaxis() -> GetBinCenter(1) - deltaTMap -> GetYaxis() -> GetBinCenter(1)) / dNegThetaZoom);

	//  Setting up waveforms and correlations to be used in interferometric map.
	TGraph waveformAnt1 = TGraph(* filtEvent -> getFilteredGraph(ant1, pol) -> even());
	TGraph waveformAnt2 = TGraph(* filtEvent -> getFilteredGraph(ant2, pol) -> even());

	TGraph crossCorr = getCorrGraph(& waveformAnt1, & waveformAnt2);
	crossCorr.SetBit(TGraph::kIsSortedX);  //  This should significantly expedite interpolation.
	TSpline3 crossCorrSpline("crossCorrSpline", & crossCorr);  //  Expedite processing by evaluating the spline once, then referencing it.

	for (int phiIdx = 1; phiIdx <= NPhi; ++phiIdx) {

		double binPhi = responseMap -> GetXaxis() -> GetBinCenter(phiIdx) * DEG2RAD;

		for (int negThetaIdx = 1; negThetaIdx <= NNegTheta; ++negThetaIdx) {

			double binNegTheta = responseMap -> GetYaxis() -> GetBinCenter(negThetaIdx) * DEG2RAD;

			double deltaT = deltaTMap ? deltaTMap -> GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx) : 0;

 			double response = crossCorrSpline.Eval(deltaT);

			int responseBinIdx = responseMap -> GetBin(phiIdx, negThetaIdx);
			double responseBinContent = responseMap -> GetBinContent(responseBinIdx);
			responseMap -> SetBinContent(responseBinIdx, responseBinContent);

			if (ant1 != ant2) responseMap -> AddBinContent(responseBinIdx, response);
			else responseMap -> AddBinContent(responseBinIdx, response / 2);
		}
	}

	//  Delete file object, then close file.
//	gROOT -> cd();
	delete deltaTMap;
//
//	gROOT -> cd();
//	deltaTMapFile.Close();
}


/*
 * Function which iteratively calls "fillMapsFlatPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillFlatMapsAll(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<vector<int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	responseMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillFlatMapsPair(filtEvent, & mapPair[i], antennaPairs[i][0], antennaPairs[i][1], pol);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) responseMap -> Add(& mapPair[i]);

	responseMap -> GetXaxis() -> SetTitle("#phi");
	responseMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * Same principal as fillCoverageMapsPair(), but no curvature included.
 * But rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillFlatCoverageMapsPair(TH2D * coverageMap, vector<vector<double>> totalPowers, int ant1, int ant2) {

	//  Access coverage maps to use.
	double totalPower1 = totalPowers[0][ant1] + totalPowers[1][ant1];
	double totalPower2 = totalPowers[0][ant2] + totalPowers[1][ant2];

	// Value to fill map with.
	double rho = sqrt(totalPower1 * totalPower2);

	for (int phiIdx = 1; phiIdx <= coverageMap -> GetNbinsX(); ++phiIdx) {

		for (int negThetaIdx = 1; negThetaIdx <= coverageMap -> GetNbinsY(); ++negThetaIdx) {

			int coverageBinIdx = coverageMap -> GetBin(phiIdx, negThetaIdx);
			double coverageBinContent = coverageMap -> GetBinContent(coverageBinIdx);
			coverageMap -> SetBinContent(coverageBinIdx, coverageBinContent);

			if (ant1 != ant2) coverageMap -> AddBinContent(coverageBinIdx, rho);
			else coverageMap -> AddBinContent(coverageBinIdx, rho / 2);
		}
	}
}


/*
 * Function which iteratively calls "fillFlatCoverageMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 * Rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillFlatCoverageMapsAll(TH2D * coverageMap, vector<vector<double>> totalPowers, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<vector<int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair coverage.
	vector<TH2D> mapPair(antennaPairs.size());
	coverageMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair coverage histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillFlatCoverageMapsPair(& mapPair[i], totalPowers, antennaPairs[i][0], antennaPairs[i][1]);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) coverageMap -> Add(& mapPair[i]);

	coverageMap -> GetXaxis() -> SetTitle("#phi");
	coverageMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 * No normalization applied. Unlike makePeakUnnormalizedInterferometricMap(), uses flat maps.
 */
TH2D makePeakUnnormalizedFlatInterferometricMap(TH2D * responseMap, int eventNum, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  First, determine location of the maximum peak in the interferometric map.
	vector<double> peakLocation = getInterferometricPeakLocation(responseMap);

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakLocation[0]);

	//  Set dimensions of interferometric maps.
	double minPhi = peakLocation[0] - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakLocation[0] + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakLocation[1] - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakLocation[1] + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	TString responseMapName = TString(responseMap -> GetName());
	const char * polType;
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) polType = "HPol";
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) polType = "VPol";
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakFlatResponseMap("peakFlatResponseMap", "Flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakFlatResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakFlatResponseMap;
}


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakUnnormalizedFlatInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakPhi);

	//  Set dimensions of interferometric maps.
	double minPhi = peakPhi - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakPhi + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakNegTheta - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakNegTheta + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	const char * polType;
	if (pol == AnitaPol::kHorizontal) polType = "HPol";
	else if (pol == AnitaPol::kVertical) polType = "VPol";
	else if (pol == AnitaPol::kNotAPol) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakFlatResponseMap("peakFlatResponseMap", "Flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, pol, neighboringAntennas);
	else {

		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakFlatResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakFlatResponseMap;
}


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferometric map about the largest peak in the input map.
 * Unlike makePeakInterferometricMap(), uses flat maps.
 */
TH2D makePeakFlatInterferometricMap(TH2D * responseMap, int eventNum, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "deconv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  First, determine location of the maximum peak in the interferometric map.
	vector<double> peakLocation = getInterferometricPeakLocation(responseMap);

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakLocation[0]);

	//  Set dimensions of interferometric maps.
	double minPhi = peakLocation[0] - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakLocation[0] + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakLocation[1] - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakLocation[1] + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	TString responseMapName = TString(responseMap -> GetName());
	const char * polType;
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) polType = "HPol";
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) polType = "VPol";
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakFlatResponseMap("peakFlatResponseMap", "Flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatCoverageMap("peakFlatCoverageMap", "Antenna pair flat coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatNormalizedResponseMap("peakFlatNormalizedResponseMap", TString::Format("%s flat normalized response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<vector<double>> totalPowers = getEventTotalPowers(& filtEvent);

	fillFlatCoverageMapsAll(& peakFlatCoverageMap, totalPowers, neighboringAntennas);

	peakFlatNormalizedResponseMap.Add(& peakFlatResponseMap);
	peakFlatNormalizedResponseMap.Divide(& peakFlatCoverageMap);
	if (!responseMapName.Contains("Unpol", TString::kIgnoreCase)) peakFlatNormalizedResponseMap.Scale(2);

	peakFlatNormalizedResponseMap.GetXaxis() -> SetTitle("#phi");
	peakFlatNormalizedResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakFlatNormalizedResponseMap;
}


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakFlatInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, TString filterString, int anitaVer, int iceMCRun) {

	AnitaVersion::set(anitaVer);

	FilterStrategy fStrat;
	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::BH13Filter());	//  Handling problematic A4 antennas.
//	UCorrelator::fillStrategyWithKey(& fStrat, "decohttps://www.netflix.com/browsenv");
//	if (anitaVer == 4) fStrat.addOperation(new UCorrelator::AntiBH13Filter());
	UCorrelator::fillStrategyWithKey(& fStrat, filterString.Data());  //  Sine subtraction filter strategy used throughout UCorrelator.

	//  Further setting up of analyzer given event number.
	int run = iceMCRun ? iceMCRun : AnitaDataset::getRunContainingEventNumber(eventNum);

	//  For when using iceMC.
	AnitaDataset::DataDirectory dataDirectory = iceMCRun ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA;

	AnitaDataset d(run, false, WaveCalType::kDefault, dataDirectory, AnitaDataset::kRandomizePolarity);
	d.getEvent(eventNum);  //  In the run, pointing to the correct eventNumber.
	FilteredAnitaEvent filtEvent(d.useful(), & fStrat, d.gps(), d.header());

	gStyle -> SetOptStat(0);  //  To remove the legend reporting number of bins.

	//  Next, determine neighboring antennas around maximum peak.
	vector<int> neighboringAntennas = getNeighboringAntennas(peakPhi);

	//  Set dimensions of interferometric maps.
	double minPhi = peakPhi - NPhiZoom * dPhiZoom / 2;
	double maxPhi = peakPhi + NPhiZoom * dPhiZoom / 2;
	double minNegTheta = peakNegTheta - NNegThetaZoom * dNegThetaZoom / 2;
	double maxNegTheta = peakNegTheta + NNegThetaZoom * dNegThetaZoom / 2;

	//  Check name of input interferometric map.
	const char * polType;
	if (pol == AnitaPol::kHorizontal) polType = "HPol";
	else if (pol == AnitaPol::kVertical) polType = "VPol";
	else if (pol == AnitaPol::kNotAPol) polType = "Unpolarized";

	//  Create TH2D objects from which to generate the peak interferometric map.
	TH2D peakFlatResponseMap("peakFlatResponseMap", "Flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatCoverageMap("peakFlatCoverageMap", "Antenna pair flat coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatNormalizedResponseMap("peakFlatNormalizedResponseMap", TString::Format("%s flat normalized response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, pol, neighboringAntennas);
	else {

		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillFlatMapsAll(& filtEvent, & peakFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<vector<double>> totalPowers = getEventTotalPowers(& filtEvent);

	fillFlatCoverageMapsAll(& peakFlatCoverageMap, totalPowers, neighboringAntennas);

	peakFlatNormalizedResponseMap.Add(& peakFlatResponseMap);
	peakFlatNormalizedResponseMap.Divide(& peakFlatCoverageMap);
	if (pol != AnitaPol::kNotAPol) peakFlatNormalizedResponseMap.Scale(2);

	peakFlatNormalizedResponseMap.GetXaxis() -> SetTitle("#phi");
	peakFlatNormalizedResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakFlatNormalizedResponseMap;
}
