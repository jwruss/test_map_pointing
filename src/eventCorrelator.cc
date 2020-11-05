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


/*
 * Return coarsely-binned histogram of time delays between antenna baselines.
 */
TH2F * getDeltaTCoarse(const char * mapName) {

	static std::map<TString, TH2F *> deltaTCoarseMap;

	if (deltaTCoarseMap.empty()) {

		#pragma omp critical
		if (deltaTCoarseMap.empty()) {  //  According to Cosmin, the purpose of the double guard is to first check if another thread is filling the map, then to ensure that it has been filled in the current thread.

			TFile deltaTCoarseFile("deltaTCoarse.root");
			TIter keyList(deltaTCoarseFile.GetListOfKeys());    //  Method of looping over objects in ROOT binary from https://root.cern.ch/root/htmldoc/tutorials/io/loopdir.C.html
			TKey * key;

			while ((key = (TKey *) keyList())) {

				TH2F * deltaT = (TH2F *) key -> ReadObj();
				deltaT -> SetDirectory(0);

				deltaTCoarseMap.insert(make_pair(key -> GetName(), deltaT));  //  This method of loading filling in the map from https://thispointer.com/stdmap-tutorial-part-1-usage-detail-with-examples/
			}
		}
	}

	return deltaTCoarseMap[mapName];
}


/*
 * Return finely-binned histogram of time delays between antenna baselines.
 */
TH2F * getDeltaTFine(const char * mapName) {

	static std::map<TString, TH2F *> deltaTFineMap;

	if (deltaTFineMap.empty()) {

		#pragma omp critical
		if (deltaTFineMap.empty()) {

			TFile deltaTFineFile("deltaTFine.root");
			TIter keyList(deltaTFineFile.GetListOfKeys());
			TKey * key;

			while ((key = (TKey *) keyList())) {

				TH2F * deltaT = (TH2F *) key -> ReadObj();
				deltaT -> SetDirectory(0);

				deltaTFineMap.insert(make_pair(key -> GetName(), deltaT));
			}
		}
	}

	return deltaTFineMap[mapName];
}


/*
 * Return coarsely-binned histogram of spherical cosine products between antenna baselines.
 */
TH2F * getSphCosProductCoarse(const char * mapName) {

	static std::map<TString, TH2F *> sphCosProductCoarseMap;

	if (sphCosProductCoarseMap.empty()) {

		#pragma omp critical
		if (sphCosProductCoarseMap.empty()) {

			TFile sphCosProductCoarseFile("sphCosProductCoarse.root");
			TIter keyList(sphCosProductCoarseFile.GetListOfKeys());
			TKey * key;

			while ((key = (TKey *) keyList())) {

				TH2F * sphCosProduct = (TH2F *) key -> ReadObj();
				sphCosProduct -> SetDirectory(0);

				sphCosProductCoarseMap.insert(make_pair(key -> GetName(), sphCosProduct));
			}
		}
	}

	return sphCosProductCoarseMap[mapName];
}


/*
 * Return finely-binned histogram of spherical cosine products between antenna baselines.
 */
TH2F * getSphCosProductFine(const char * mapName) {

	static std::map<TString, TH2F *> sphCosProductFineMap;

	if (sphCosProductFineMap.empty()) {

		#pragma omp critical
		if (sphCosProductFineMap.empty()) {

			TFile sphCosProductFineFile("sphCosProductFine.root");
			TIter keyList(sphCosProductFineFile.GetListOfKeys());
			TKey * key;

			while ((key = (TKey *) keyList())) {

				TH2F * sphCosProduct = (TH2F *) key -> ReadObj();
				sphCosProduct -> SetDirectory(0);

				sphCosProductFineMap.insert(make_pair(key -> GetName(), sphCosProduct));
			}
		}
	}

	return sphCosProductFineMap[mapName];
}


/*
 * Get antenna pairs used to construct interferometric maps.
 */
vector<pair<int, int>> getAntennaPairs(vector<int> neighboringAntennas) {

	//  Check if all neighboring antennas allowed.
	if (neighboringAntennas.empty()) {

		neighboringAntennas.resize(NUM_SEAVEYS);
		for (int i = 0; i < NUM_SEAVEYS; ++i) neighboringAntennas[i] = i;
	}

	vector<pair<int, int>> antennaPairs;  //  Where antenna pairs will be placed.

	//  Nested for loop which fills the antennaPairs vector.
	for (int i = 0; i < neighboringAntennas.size(); ++i) {

		for (int j = i + 1; j < neighboringAntennas.size(); ++j) {

			int phiSep = abs(neighboringAntennas[i] - neighboringAntennas[j]) % NUM_PHI;
			phiSep = min(phiSep, NUM_PHI - phiSep);
			if (phiSep > 2) continue;  //  Exclude pairs more than 2 phi sectors apart.

			antennaPairs.push_back(make_pair(neighboringAntennas[i], neighboringAntennas[j]));
		}
	}

	return antennaPairs;
}


/*
 * Given an event, construct the cross-correlation between two antennas for a given polarization.
 */
void fillPairMap(FilteredAnitaEvent * filtEvent, TH2D * responseMap, pair<int, int> antPair, AnitaPol::AnitaPol_t pol) {

	//  Check dimensions of responseMap, to determine which spherical cosine pair file to reference.
	int NPhi = responseMap -> GetNbinsX();
	int NNegTheta = responseMap -> GetNbinsY();

	//  Char associated with polarization.
	const char * polChar = (pol == AnitaPol::kHorizontal) ? "H" : "V";

	//  Determine which pair of spherical cosines histogram to reference.
	TH2F * deltaTMap = 0;
	if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

		deltaTMap = getDeltaTCoarse(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());
//		deltaTMap -> Scale(-1);  //  Testing reciprocity. Normally this is commented out.

		if (!deltaTMap) {  //  Check if NULL pointer.

			deltaTMap = getDeltaTCoarse(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());
			deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
		}

	} else {

		deltaTMap = getDeltaTFine(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());

		if (!deltaTMap) {

			deltaTMap = getDeltaTFine(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());
			deltaTMap -> Scale(-1);
		}
	}

//	deltaTMap -> SetDirectory(0);

	TH2F * sphCosProduct = 0;
	if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

		sphCosProduct = getSphCosProductCoarse(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());

		if (!sphCosProduct) sphCosProduct= getSphCosProductCoarse(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());

	} else {

		sphCosProduct = getSphCosProductFine(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());

		if (!sphCosProduct) sphCosProduct = getSphCosProductFine(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());
	}

//	sphCosProduct -> SetDirectory(0);

//	//  Accessing relevant TH2F objects.
//	TFile deltaTFile((NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? "deltaTCoarse.root" : "deltaTFine.root");
//	TH2F * deltaTMap = (TH2F *) deltaTFile.Get(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second));
//	if (!deltaTMap) {
//
//		deltaTMap = (TH2F *) deltaTFile.Get(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first));
//		deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
//	}
////	deltaTMap -> SetDirectory(0);
//
//	TFile sphCosProductFile((NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? "sphCosProductCoarse.root" : "sphCosProductFine.root");
//	TH2F * sphCosProduct = (TH2F *) sphCosProductFile.Get(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second));
//	if (!sphCosProduct) sphCosProduct = (TH2F *) sphCosProductFile.Get(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first));
////	sphCosProduct -> SetDirectory(0);

	//  Get index differences between reseponseMap and sphCosProduct. Relevant for zoomed maps.
	int dPhiIdx = int((responseMap -> GetXaxis() -> GetBinCenter(1) - deltaTMap -> GetXaxis() -> GetBinCenter(1)) / dPhiZoom);
	int dNegThetaIdx = int((responseMap -> GetYaxis() -> GetBinCenter(1) - deltaTMap -> GetYaxis() -> GetBinCenter(1)) / dNegThetaZoom);

	//  Setting up waveforms and correlations to be used in interferometric map.
	TGraph waveformAnt1 = TGraph(* filtEvent -> getFilteredGraph(antPair.first, pol) -> even());
	TGraph waveformAnt2 = TGraph(* filtEvent -> getFilteredGraph(antPair.second, pol) -> even());

	TGraph crossCorr = getCorrGraph(& waveformAnt1, & waveformAnt2);
	crossCorr.SetBit(TGraph::kIsSortedX);  //  This should significantly expedite interpolation.
	TSpline3 crossCorrSpline("crossCorrSpline", & crossCorr);  //  Expedite processing by evaluating the spline once, then referencing it.

	for (int phiIdx = 1; phiIdx <= NPhi; ++phiIdx) {

		double binPhi = responseMap -> GetXaxis() -> GetBinCenter(phiIdx) * DEG2RAD;

		for (int negThetaIdx = 1; negThetaIdx <= NNegTheta; ++negThetaIdx) {

			double rho = sphCosProduct -> GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx);

			if (rho <= 0) continue;  //  Reference the spherical cosine pair to determine if bin should be filled.

			double binNegTheta = responseMap -> GetYaxis() -> GetBinCenter(negThetaIdx) * DEG2RAD;

			double deltaT = deltaTMap ? deltaTMap -> GetBinContent(phiIdx + dPhiIdx, negThetaIdx + dNegThetaIdx) : 0;

 			double response = crossCorrSpline.Eval(deltaT);

			int responseBinIdx = responseMap -> GetBin(phiIdx, negThetaIdx);
			double responseBinContent = responseMap -> GetBinContent(responseBinIdx);
			responseMap -> SetBinContent(responseBinIdx, responseBinContent);

			if (antPair.first != antPair.second) responseMap -> AddBinContent(responseBinIdx, rho * response);
			else responseMap -> AddBinContent(responseBinIdx, rho * response / 2);
		}
	}

//	//  Delete pointers to file objects, then close the files.
//	gROOT -> cd();
//	delete deltaTMap;
//	delete sphCosProduct;
//
//	gROOT -> cd();
//	deltaTMapFile.Close();
//	antSphericalCosineProductsFile.Close();
}


/*
 * Function which iteratively calls "fillMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllMaps(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<pair<int, int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	responseMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillPairMap(filtEvent, & mapPair[i], antennaPairs[i], pol);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) responseMap -> Add(& mapPair[i]);

	responseMap -> GetXaxis() -> SetTitle("#phi");
	responseMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * For a given event, return a 2-dimensional vector containing sums of square voltages for each antenna. Vector has dimensions of
 * [2][NUM_SEAVEYS]. where first index corresponds to 0 for Hpol, 1 for Vpol.
 */
vector<pair<double, double>> getEventTotalPowers(FilteredAnitaEvent * filtEvent) {

	vector<pair<double, double>> totalPower(NUM_SEAVEYS);

	#pragma omp parallel for
	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		TGraph waveformAntH = TGraph(* filtEvent -> getFilteredGraph(i, AnitaPol::kHorizontal) -> even());
		TGraph waveformAntV = TGraph(* filtEvent -> getFilteredGraph(i, AnitaPol::kVertical) -> even());

		for (int j = 0; j < waveformAntH.GetN(); ++j) totalPower[i].first += waveformAntH.GetY()[j] * waveformAntH.GetY()[j];
		for (int j = 0; j < waveformAntV.GetN(); ++j) totalPower[i].second += waveformAntV.GetY()[j] * waveformAntV.GetY()[j];
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
	vector<pair<double, double>> totalPowers = getEventTotalPowers(filtEvent);

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
//	TFile sphCosProductFile(isFine ? "sphCosProductFine.root" : "sphCosProductCoarse.root");

	#pragma omp parallel for
	for (int i = 0; i < NUM_SEAVEYS; ++i) {

		TH2F * hHist = 0, * vHist = 0;

		if (isFine) {

			hHist = getSphCosProductFine(TString::Format("H_%d_%d", i, i).Data());
			vHist = getSphCosProductFine(TString::Format("V_%d_%d", i, i).Data());

		} else {

			hHist = getSphCosProductCoarse(TString::Format("H_%d_%d", i, i).Data());
			vHist = getSphCosProductCoarse(TString::Format("V_%d_%d", i, i).Data());
		}

//		hHist -> SetDirectory(0);
//		vHist -> SetDirectory(0);

//		TH2F * hHist = (TH2F *) sphCosProductFile.Get(TString::Format("H_%d_%d", i, i));
////		hHist -> SetDirectory(0);
//
//		TH2D * vHist = (TH2F *) sphCosProductFile.Get(TString::Format("V_%d_%d", i, i));
////		vHist -> SetDirectory(0);

		totalPowerMaps[i] = TH2D(TString::Format("Ant_%d", i), TString::Format("Ant_%d", i), NPhi, minPhi, maxPhi, NNegTheta, minNegTheta, maxNegTheta);
		totalPowerMaps[i].Add(hHist, vHist, totalPowers[i].first, totalPowers[i].first);

		totalPowerMaps[i].GetXaxis() -> SetTitle("#phi");
		totalPowerMaps[i].GetYaxis() -> SetTitle("-#theta");

//		gROOT -> cd();
//		delete hHist;
//		delete vHist;
	}

//	//  Close the file.
//	gROOT -> cd();
//	antSphericalCosineProductsFile.Close();

	return totalPowerMaps;
}


/*
 * Each antenna cross-correlation has an upper bound given by the Cauchy-Schwarz inequality. The upper bound here is over both horizontal and vertical polarization.
 */
void fillPairCoverageMap(TH2D * coverageMap, vector<TH2D> totalPowerMaps, pair<int, int> antPair) {

	//  Access coverage maps to use.
	TH2D totalPowerMap1 = totalPowerMaps[antPair.first];
	TH2D totalPowerMap2 = totalPowerMaps[antPair.second];

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

			if (antPair.first != antPair.second) coverageMap -> AddBinContent(coverageBinIdx, rho);
			else coverageMap -> AddBinContent(coverageBinIdx, rho / 2);
		}
	}
}


/*
 * Function which iteratively calls "fillCoverageMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllCoverageMaps(TH2D * coverageMap, vector<TH2D> totalPowerMaps, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<pair<int, int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	coverageMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair coverage histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillPairCoverageMap(& mapPair[i], totalPowerMaps, antennaPairs[i]);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) coverageMap -> Add(& mapPair[i]);

	coverageMap -> GetXaxis() -> SetTitle("#phi");
	coverageMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * For a given event, calls to "fillAllMaps()" to construct corresponding unnormalized interferoemtric maps.
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
	fillAllMaps(& filtEvent, & responseMapHPol, AnitaPol::kHorizontal);
	fillAllMaps(& filtEvent, & responseMapVPol, AnitaPol::kVertical);

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
 * For a given event, calls to "fillAllMaps()" and "fillCoverageMapsAll()" to construct corresponding interferometric maps.
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
	fillAllMaps(& filtEvent, & responseMapHPol, AnitaPol::kHorizontal);
	fillAllMaps(& filtEvent, & responseMapVPol, AnitaPol::kVertical);

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent);

	fillAllCoverageMaps(& coverageMap, totalPowerMaps);

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
	TH2D peakUnnormalizedResponseMap(!useBroadband ? "peakUnnormalizedResponseMap" : "peakUnnormalizedBroadbandResponseMap", !useBroadband ? "Unnormalized response around interferometric peak" : "Unnormalized broadband response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakUnnormalizedResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakUnnormalizedResponseMap;
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
	TH2D peakUnnormalizedResponseMap(!useBroadband ? "peakUnnormalizedResponseMap" : "peakUnnormalizedBroadbandResponseMap", !useBroadband ? "Unnormalized response around interferometric peak" : "Unnormalized broadband response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, pol, neighboringAntennas);
	else {

		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakUnnormalizedResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakUnnormalizedResponseMap;
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
	TH2D peakUnnormalizedResponseMap("peakUnnormalizedResponseMap", "Unnormalized response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakCoverageMap("peakCoverageMap", "Antenna pair coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakResponseMap(!useBroadband ? "peakResponseMap" : "peakBroadbandResponseMap", !useBroadband ? TString::Format("%s response around interferometric peak", polType) : TString::Format("%s broadband response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillAllMaps(& filtEvent, & peakResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillAllMaps(& filtEvent, & peakResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent, true);

	fillAllCoverageMaps(& peakCoverageMap, totalPowerMaps, neighboringAntennas);

	peakResponseMap.Add(& peakUnnormalizedResponseMap);
	peakResponseMap.Divide(& peakCoverageMap);
	if (!responseMapName.Contains("Unpol", TString::kIgnoreCase)) peakResponseMap.Scale(2);

	peakResponseMap.GetXaxis() -> SetTitle("#phi");
	peakResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakResponseMap;
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
	TH2D peakUnnormalizedResponseMap("peakUnnormalizedResponseMap", "Unnormalized response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakCoverageMap("peakCoverageMap", "Antenna pair coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakResponseMap(!useBroadband ? "peakResponseMap" : "peakBroadbandResponseMap", !useBroadband ? TString::Format("%s response around interferometric peak", polType) : TString::Format("%s broadband response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, pol, neighboringAntennas);
	else {

		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllMaps(& filtEvent, & peakUnnormalizedResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<TH2D> totalPowerMaps = getTotalPowerMaps(& filtEvent, true);

	fillAllCoverageMaps(& peakCoverageMap, totalPowerMaps, neighboringAntennas);

	peakResponseMap.Add(& peakUnnormalizedResponseMap);
	peakResponseMap.Divide(& peakCoverageMap);
	if (pol != AnitaPol::kNotAPol) peakResponseMap.Scale(2);

	peakResponseMap.GetXaxis() -> SetTitle("#phi");
	peakResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakResponseMap;
}


/*
 * For purposes of testing noncoplanar baselines, create cross-correlation
 * interferometric maps without curvature correction.
 */
void fillPairFlatMap(FilteredAnitaEvent * filtEvent, TH2D * responseMap, pair<int, int> antPair, AnitaPol::AnitaPol_t pol) {

	//  Check dimensions of responseMap, to determine which spherical cosine pair file to reference.
	int NPhi = responseMap -> GetNbinsX();
	int NNegTheta = responseMap -> GetNbinsY();

	//  Char associated with polarization.
	const char * polChar = (pol == AnitaPol::kHorizontal) ? "H" : "V";

	//  Determine which pair of spherical cosines histogram to reference.
	TH2F * deltaTMap = 0;
	if (NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) {

		deltaTMap = getDeltaTCoarse(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());

		if (!deltaTMap) {  //  Check if NULL pointer.

			deltaTMap = getDeltaTCoarse(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());
			deltaTMap -> Scale(-1);  //  Ordering is antisymmetric here.
		}

	} else {

		deltaTMap = getDeltaTFine(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second).Data());

		if (!deltaTMap) {

			deltaTMap = getDeltaTFine(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first).Data());
			deltaTMap -> Scale(-1);
		}
	}

//	deltaTMap -> SetDirectory(0);

//	TFile deltaTFile((NPhi == NPhiCoarse && NNegTheta == NNegThetaCoarse) ? "deltaTCoarse.root" : "deltaTFine.root");
//	TH2F * deltaTMap = (TH2F *) deltaTFile.Get(TString::Format("%s_%d_%d", polChar, antPair.first, antPair.second));
//	if (!deltaTMap) {
//
//		deltaTMap = (TH2F *) deltaTFile.Get(TString::Format("%s_%d_%d", polChar, antPair.second, antPair.first));
//		deltaTMap -> Scale(-1);  //  Ordering is antisymmetric, here.
//	}
////	deltaTMap -> SetDirectory(0);

	//  Get index differences between reseponseMap and sphCosProduct. Relevant for zoomed maps.
	int dPhiIdx = int((responseMap -> GetXaxis() -> GetBinCenter(1) - deltaTMap -> GetXaxis() -> GetBinCenter(1)) / dPhiZoom);
	int dNegThetaIdx = int((responseMap -> GetYaxis() -> GetBinCenter(1) - deltaTMap -> GetYaxis() -> GetBinCenter(1)) / dNegThetaZoom);

	//  Setting up waveforms and correlations to be used in interferometric map.
	TGraph waveformAnt1 = TGraph(* filtEvent -> getFilteredGraph(antPair.first, pol) -> even());
	TGraph waveformAnt2 = TGraph(* filtEvent -> getFilteredGraph(antPair.second, pol) -> even());

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

			if (antPair.first != antPair.second) responseMap -> AddBinContent(responseBinIdx, response);
			else responseMap -> AddBinContent(responseBinIdx, response / 2);
		}
	}

	//  Delete file object, then close file.
//	gROOT -> cd();
//	delete deltaTMap;
//
//	gROOT -> cd();
//	deltaTMapFile.Close();
}


/*
 * Function which iteratively calls "fillMapsFlatPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllFlatMaps(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<pair<int, int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair response.
	vector<TH2D> mapPair(antennaPairs.size());
	responseMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillPairFlatMap(filtEvent, & mapPair[i], antennaPairs[i], pol);

	//  Add the histograms together.
	for (int i = 0; i < antennaPairs.size(); ++i) responseMap -> Add(& mapPair[i]);

	responseMap -> GetXaxis() -> SetTitle("#phi");
	responseMap -> GetYaxis() -> SetTitle("-#theta");
}


/*
 * Same principal as fillCoverageMapsPair(), but no curvature included.
 * But rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillPairFlatCoverageMap(TH2D * coverageMap, vector<pair<double, double>> totalPowers, pair<int, int> antPair) {

	//  Access coverage maps to use.
	double totalPower1 = totalPowers[antPair.first].first + totalPowers[antPair.first].second;
	double totalPower2 = totalPowers[antPair.second].first + totalPowers[antPair.second].second;

	// Value to fill map with.
	double rho = sqrt(totalPower1 * totalPower2);

	for (int phiIdx = 1; phiIdx <= coverageMap -> GetNbinsX(); ++phiIdx) {

		for (int negThetaIdx = 1; negThetaIdx <= coverageMap -> GetNbinsY(); ++negThetaIdx) {

			int coverageBinIdx = coverageMap -> GetBin(phiIdx, negThetaIdx);
			double coverageBinContent = coverageMap -> GetBinContent(coverageBinIdx);
			coverageMap -> SetBinContent(coverageBinIdx, coverageBinContent);

			if (antPair.first != antPair.second) coverageMap -> AddBinContent(coverageBinIdx, rho);
			else coverageMap -> AddBinContent(coverageBinIdx, rho / 2);
		}
	}
}


/*
 * Function which iteratively calls "fillFlatCoverageMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 * Rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillAllFlatCoverageMaps(TH2D * coverageMap, vector<pair<double, double>> totalPowers, vector<int> neighboringAntennas) {

	//  Establish antenna pairs to be use.
	vector<pair<int, int>> antennaPairs = getAntennaPairs(neighboringAntennas);

	//  Create vector array in which to place each antenna pair coverage.
	vector<TH2D> mapPair(antennaPairs.size());
	coverageMap -> Copy(mapPair[0]);
	mapPair[0].Reset();
	for (int i = 1; i < antennaPairs.size(); ++i) mapPair[0].Copy(mapPair[i]);

	//  Fill the antenna pair coverage histograms. This loop should be embarrassingly parallel.
	#pragma omp parallel for
	for (int i = 0; i < antennaPairs.size(); ++i) fillPairFlatCoverageMap(& mapPair[i], totalPowers, antennaPairs[i]);

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
	TH2D peakUnnormalizedFlatResponseMap("peakUnnormalizedFlatResponseMap", "Unnormalized flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakUnnormalizedFlatResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakUnnormalizedFlatResponseMap;
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
	TH2D peakUnnormalizedFlatResponseMap("peakUnnormalizedFlatResponseMap", "Unnormalized flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, pol, neighboringAntennas);
	else {

		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
		peakUnnormalizedFlatResponseMap.Scale(0.5);
	}

	//  Return pointer to peak interferometric map.
	return peakUnnormalizedFlatResponseMap;
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
	TH2D peakUnnormalizedFlatResponseMap("peakUnnormalizedFlatResponseMap", "Unnormalized flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatCoverageMap("peakFlatCoverageMap", "Antenna pair flat coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatResponseMap("peakFlatResponseMap", TString::Format("%s flat response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (responseMapName.Contains("HPol", TString::kIgnoreCase)) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
	else if (responseMapName.Contains("VPol", TString::kIgnoreCase)) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	else if (responseMapName.Contains("Unpol", TString::kIgnoreCase)) {

		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<pair<double, double>> totalPowers = getEventTotalPowers(& filtEvent);

	fillAllFlatCoverageMaps(& peakFlatCoverageMap, totalPowers, neighboringAntennas);

	peakFlatResponseMap.Add(& peakUnnormalizedFlatResponseMap);
	peakFlatResponseMap.Divide(& peakFlatCoverageMap);
	if (!responseMapName.Contains("Unpol", TString::kIgnoreCase)) peakFlatResponseMap.Scale(2);

	peakFlatResponseMap.GetXaxis() -> SetTitle("#phi");
	peakFlatResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakFlatResponseMap;
}


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakFlatInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, TString filterString, int anitaVer, int iceMCRun) {

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
	TH2D peakUnnormalizedFlatResponseMap("peakUnnormalizedFlatResponseMap", "Unnormalized flat response around interferometric peak", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatCoverageMap("peakFlatCoverageMap", "Antenna pair flat coverage around interferometric peak.", NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);
	TH2D peakFlatResponseMap("peakFlatResponseMap", TString::Format("%s flat response around interferometric peak", polType), NPhiZoom, minPhi, maxPhi, NNegThetaZoom, minNegTheta, maxNegTheta);

	//  Use preceding functions about peak interferometric map.
	if (pol != AnitaPol::kNotAPol) fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, pol, neighboringAntennas);
	else {

		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kHorizontal, neighboringAntennas);
		fillAllFlatMaps(& filtEvent, & peakUnnormalizedFlatResponseMap, AnitaPol::kVertical, neighboringAntennas);
	}

	vector<pair<double, double>> totalPowers = getEventTotalPowers(& filtEvent);

	fillAllFlatCoverageMaps(& peakFlatCoverageMap, totalPowers, neighboringAntennas);

	peakFlatResponseMap.Add(& peakUnnormalizedFlatResponseMap);
	peakFlatResponseMap.Divide(& peakFlatCoverageMap);
	if (pol != AnitaPol::kNotAPol) peakFlatResponseMap.Scale(2);

	peakFlatResponseMap.GetXaxis() -> SetTitle("#phi");
	peakFlatResponseMap.GetYaxis() -> SetTitle("-#theta");

	//  Return pointer to peak interferometric map.
	return peakFlatResponseMap;
}
