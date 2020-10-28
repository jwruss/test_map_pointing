/*
 * testMapPointing.cxx
 *
 *  Created on: Oct 08, 2020
 *      Author: John Russell
 *
 *  With "eventCorrelator.h", here we construct pointing in recalculated,
 *  zoomed interferometric maps. Using WAIS events, this is to confirm whether
 *  or not the changes I would like to make to how we calculate interferometric
 *  maps indeed helps with pointing resolution.
 */

using namespace std;

#include "AnitaDataset.h"
#include "AnitaConventions.h"
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"
#include "AnitaEventSummary.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "eventCorrelator.h"
//#include "TROOT.h"

const double dPhiRough = 2.;
const double dThetaRough = 1.;


void addPolTree(int part, bool isVPol = false) {

	//  First declare a const char * associated with chose polarization.
	const char * polChar = !isVPol ? "H" : "V";

	//  Event summary file from which to reference events for given run.
	TFile sumFile(TString::Format("~/jwrussWork/eventSummaries/WAIS%sPol/summary-full-WAIS%sPol-part_%d.root", polChar, polChar, part));

	//  Point to relevant quantities in the summary file.
	TTree * sumTree = (TTree *) sumFile.Get("sampleA4");
	AnitaEventSummary * sum = 0;
	sumTree -> SetBranchAddress("summary", & sum);

	//  Reopen file in which to add TTree.
	TFile pointingFile(TString::Format("testMapPointingFiles/testMapPointing-part_%d.root", part), "update");

	//  Creating TTree in which to store pointing results.
	TTree pointingTree(TString::Format("WAIS%sPolPointingTree", polChar), TString::Format("TTree containing pointing information for %s polarized WAIS events.", !isVPol ? "horizontally" : "vertically"));
	int run;
	pointingTree.Branch("run", & run);
	unsigned int eventNumber;
	pointingTree.Branch("eventNumber", & eventNumber);
	double cohLinPolAngle, deconvLinPolAngle;
	pointingTree.Branch("cohLinPolAngle", & cohLinPolAngle);
	pointingTree.Branch("deconvLinPolAngle", & deconvLinPolAngle);
	double phiRough, thetaRough, phiWAIS, thetaWAIS;  //  Initial values coming from AnitaEventSummary.
	pointingTree.Branch("phiRough", & phiRough);
	pointingTree.Branch("thetaRough", & thetaRough);
	pointingTree.Branch("phiWAIS", & phiWAIS);
	pointingTree.Branch("thetaWAIS", & thetaWAIS);
	double peakDefault, phiDefault, thetaDefault;  //  The first set are initial v
	pointingTree.Branch("peakDefault", & peakDefault);
	pointingTree.Branch("phiDefault", & phiDefault);
	pointingTree.Branch("thetaDefault", & thetaDefault);

	struct {

		double peakPol, phiPol, thetaPol;
		double peakUnpol, phiUnpol, thetaUnpol;

	} unnormalizedFlat, unnormalized, unnormalizedBroadband, normalizedFlat, normalized, normalizedBroadband;

	TString structStr = "peakPol/D:phiPol:thetaPol:peakUnpol:phiUnpol:thetaUnpol";
	pointingTree.Branch("unnormalizedFlat", & unnormalizedFlat, structStr);
	pointingTree.Branch("unnormalized", & unnormalized, structStr);
	pointingTree.Branch("unnormalizedBroadband", & unnormalizedBroadband, structStr);
	pointingTree.Branch("normalizedFlat", & normalizedFlat, structStr);
	pointingTree.Branch("normalized", & normalized, structStr);
	pointingTree.Branch("normalizedBroadband", & normalizedBroadband, structStr);

	//  Iterate through the summary file, creating finely binned interferometric maps and extracting peak location and value.
//	int numEntries = sumTree -> GetEntries() > 15 ? 15 : sumTree -> GetEntries();
	int numEntries = sumTree -> GetEntries();
	for (int entryNum = 0; entryNum < numEntries; ++entryNum) {

		sumTree -> GetEntry(entryNum);

		run = sum -> run;
		eventNumber = sum -> eventNumber;

		cohLinPolAngle = sum -> highestCoherentFiltered().linearPolAngle();
		deconvLinPolAngle = sum -> highestDeconvolvedFiltered().linearPolAngle();
//		cohLinPolAngle = sum -> mostImpulsiveCoherentFiltered(2).linearPolAngle();
//		deconvLinPolAngle = sum -> mostImpulsiveDeconvolvedFiltered(2).linearPolAngle();

		phiRough = dPhiRough * (ceil(FFTtools::wrap(sum -> highestPeak().phi + sum -> highestPeak().dphi_rough) / dPhiRough) - 0.5);
		thetaRough = dThetaRough * (ceil((sum -> highestPeak().theta + sum -> highestPeak().dtheta_rough) / dThetaRough) - 0.5);
//		phiRough = dPhiRough * (ceil(FFTtools::wrap(sum -> mostImpulsivePeak(2).phi + sum -> mostImpulsivePeak(2).dphi_rough) / dPhiRough) - 0.5);
//		thetaRough = dThetaRough * (ceil((sum -> mostImpulsivePeak(2).theta + sum -> mostImpulsivePeak(2).dtheta_rough) / dThetaRough) - 0.5);
		phiWAIS = sum -> wais.phi;
		thetaWAIS = sum -> wais.theta;

		peakDefault = sum -> highestPeak().value;
		phiDefault = sum -> highestPeak().phi;
		thetaDefault = sum -> highestPeak().theta;
//		peakDefault = sum -> mostImpulsivePeak(2).value;
//		phiDefault = sum -> mostImpulsivePeak(2).phi;
//		thetaDefault = sum -> mostImpulsivePeak(2).theta;

//		TH2D unnormalizedFlatMapPol = TH2D(makePeakUnnormalizedFlatInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
//		unnormalizedFlat.peakPol = unnormalizedFlatMapPol.GetMaximum();
//		vector<double> unnormalizedFlatMapPolAngles = getInterferometricPeakLocation(& unnormalizedFlatMapPol);
		TH2D * unnormalizedFlatMapPol = new TH2D(makePeakUnnormalizedFlatInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
		unnormalizedFlat.peakPol = unnormalizedFlatMapPol -> GetMaximum();
		vector<double> unnormalizedFlatMapPolAngles = getInterferometricPeakLocation(unnormalizedFlatMapPol);
		unnormalizedFlat.phiPol = unnormalizedFlatMapPolAngles[0];
		unnormalizedFlat.thetaPol = -unnormalizedFlatMapPolAngles[1];
		delete unnormalizedFlatMapPol;
//		delete gROOT -> FindObject("peakFlatResponseMap");
//		TH2D unnormalizedFlatMapUnpol = TH2D(makePeakUnnormalizedFlatInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
//		unnormalizedFlat.peakUnpol = unnormalizedFlatMapUnpol.GetMaximum();
//		vector<double> unnormalizedFlatMapUnpolAngles = getInterferometricPeakLocation(& unnormalizedFlatMapUnpol);
		TH2D * unnormalizedFlatMapUnpol = new TH2D(makePeakUnnormalizedFlatInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
		unnormalizedFlat.peakUnpol = unnormalizedFlatMapUnpol -> GetMaximum();
		vector<double> unnormalizedFlatMapUnpolAngles = getInterferometricPeakLocation(unnormalizedFlatMapUnpol);
		unnormalizedFlat.phiUnpol = unnormalizedFlatMapUnpolAngles[0];
		unnormalizedFlat.thetaUnpol = -unnormalizedFlatMapUnpolAngles[1];
		delete unnormalizedFlatMapUnpol;
//		delete gROOT -> FindObject("peakFlatResponseMap");

//		TH2D unnormalizedMapPol = TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
//		unnormalized.peakPol = unnormalizedMapPol.GetMaximum();
//		vector<double> unnormalizedMapPolAngles = getInterferometricPeakLocation(& unnormalizedMapPol);
		TH2D * unnormalizedMapPol = new TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
		unnormalized.peakPol = unnormalizedMapPol -> GetMaximum();
		vector<double> unnormalizedMapPolAngles = getInterferometricPeakLocation(unnormalizedMapPol);
		unnormalized.phiPol = unnormalizedMapPolAngles[0];
		unnormalized.thetaPol = -unnormalizedMapPolAngles[1];
		delete unnormalizedMapPol;
//		delete gROOT -> FindObject("peakResponseMap");
//		TH2D unnormalizedMapUnpol = TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
//		unnormalized.peakUnpol = unnormalizedMapUnpol.GetMaximum();
//		vector<double> unnormalizedMapUnpolAngles = getInterferometricPeakLocation(& unnormalizedMapUnpol);
		TH2D * unnormalizedMapUnpol = new TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
		unnormalized.peakUnpol = unnormalizedMapUnpol -> GetMaximum();
		vector<double> unnormalizedMapUnpolAngles = getInterferometricPeakLocation(unnormalizedMapUnpol);
		unnormalized.phiUnpol = unnormalizedMapUnpolAngles[0];
		unnormalized.thetaUnpol = -unnormalizedMapUnpolAngles[1];
		delete unnormalizedMapUnpol;
//		delete gROOT -> FindObject("peakResponseMap");

//		TH2D unnormalizedBroadbandMapPol = makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> mostImpulsivePol(2), true);
//		unnormalizedBroadband.peakPol = unnormalizedBroadbandMapPol.GetMaximum();
//		vector<double> unnormalizedBroadbandMapPolAngles = getInterferometricPeakLocation(& unnormalizedBroadbandMapPol);
		TH2D * unnormalizedBroadbandMapPol = new TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol(), true));
		unnormalizedBroadband.peakPol = unnormalizedBroadbandMapPol -> GetMaximum();
		vector<double> unnormalizedBroadbandMapPolAngles = getInterferometricPeakLocation(unnormalizedBroadbandMapPol);
		unnormalizedBroadband.phiPol = unnormalizedBroadbandMapPolAngles[0];
		unnormalizedBroadband.thetaPol = -unnormalizedBroadbandMapPolAngles[1];
		delete unnormalizedBroadbandMapPol;
//		delete gROOT -> FindObject("peakResponseMap");
//		TH2D unnormalizedBroadbandMapUnpol = makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol, true);
//		unnormalizedBroadband.peakUnpol = unnormalizedBroadbandMapUnpol.GetMaximum();
//		vector<double> unnormalizedBroadbandMapUnpolAngles = getInterferometricPeakLocation(& unnormalizedBroadbandMapUnpol);
		TH2D * unnormalizedBroadbandMapUnpol = new TH2D(makePeakUnnormalizedInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol, true));
		unnormalizedBroadband.peakUnpol = unnormalizedBroadbandMapUnpol -> GetMaximum();
		vector<double> unnormalizedBroadbandMapUnpolAngles = getInterferometricPeakLocation(unnormalizedBroadbandMapUnpol);
		unnormalizedBroadband.phiUnpol = unnormalizedBroadbandMapUnpolAngles[0];
		unnormalizedBroadband.thetaUnpol = -unnormalizedBroadbandMapUnpolAngles[1];
		delete unnormalizedBroadbandMapUnpol;
//		delete gROOT -> FindObject("peakResponseMap");

//		TH2D normalizedFlatMapPol = makePeakFlatInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> mostImpulsivePol(2));
//		normalizedFlat.peakPol = normalizedFlatMapPol.GetMaximum();
//		vector<double> normalizedFlatMapPolAngles = getInterferometricPeakLocation(& normalizedFlatMapPol);
		TH2D * normalizedFlatMapPol = new TH2D(makePeakFlatInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
		normalizedFlat.peakPol = normalizedFlatMapPol -> GetMaximum();
		vector<double> normalizedFlatMapPolAngles = getInterferometricPeakLocation(normalizedFlatMapPol);
		normalizedFlat.phiPol = normalizedFlatMapPolAngles[0];
		normalizedFlat.thetaPol = -normalizedFlatMapPolAngles[1];
		delete normalizedFlatMapPol;
//		delete gROOT -> FindObject("peakFlatNormalizedResponseMap");
//		TH2D normalizedFlatMapUnpol = makePeakFlatInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol);
//		normalizedFlat.peakUnpol = normalizedFlatMapUnpol.GetMaximum();
//		vector<double> normalizedFlatMapUnpolAngles = getInterferometricPeakLocation(& normalizedFlatMapUnpol);
		TH2D * normalizedFlatMapUnpol = new TH2D(makePeakFlatInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
		normalizedFlat.peakUnpol = normalizedFlatMapUnpol -> GetMaximum();
		vector<double> normalizedFlatMapUnpolAngles = getInterferometricPeakLocation(normalizedFlatMapUnpol);
		normalizedFlat.phiUnpol = normalizedFlatMapUnpolAngles[0];
		normalizedFlat.thetaUnpol = -normalizedFlatMapUnpolAngles[1];
		delete normalizedFlatMapUnpol;
//		delete gROOT -> FindObject("peakFlatNormalizedResponseMap");

//		TH2D normalizedMapPol = makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> mostImpulsivePol(2));
//		normalized.peakPol = normalizedMapPol.GetMaximum();
//		vector<double> normalizedMapPolAngles = getInterferometricPeakLocation(& normalizedMapPol);
		TH2D * normalizedMapPol = new TH2D(makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol()));
		normalized.peakPol = normalizedMapPol -> GetMaximum();
		vector<double> normalizedMapPolAngles = getInterferometricPeakLocation(normalizedMapPol);
		normalized.phiPol = normalizedMapPolAngles[0];
		normalized.thetaPol = -normalizedMapPolAngles[1];
		delete normalizedMapPol;
//		delete gROOT -> FindObject("peakNormalizedResponseMap");
//		TH2D normalizedMapUnpol = makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol);
//		normalized.peakUnpol = normalizedMapUnpol.GetMaximum();
//		vector<double> normalizedMapUnpolAngles = getInterferometricPeakLocation(& normalizedMapUnpol);
		TH2D * normalizedMapUnpol = new TH2D(makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol));
		normalized.peakUnpol = normalizedMapUnpol -> GetMaximum();
		vector<double> normalizedMapUnpolAngles = getInterferometricPeakLocation(normalizedMapUnpol);
		normalized.phiUnpol = normalizedMapUnpolAngles[0];
		normalized.thetaUnpol = -normalizedMapUnpolAngles[1];
		delete normalizedMapUnpol;
//		delete gROOT -> FindObject("peakNormalizedResponseMap");

//		TH2D normalizedBroadbandMapPol = makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> mostImpulsivePol(2), true);
//		normalizedBroadband.peakPol = normalizedBroadbandMapPol.GetMaximum();
//		vector<double> normalizedBroadbandMapPolAngles = getInterferometricPeakLocation(& normalizedBroadbandMapPol);
		TH2D * normalizedBroadbandMapPol = new TH2D(makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, sum -> highestPol(), true));
		normalizedBroadband.peakPol = normalizedBroadbandMapPol -> GetMaximum();
		vector<double> normalizedBroadbandMapPolAngles = getInterferometricPeakLocation(normalizedBroadbandMapPol);
		normalizedBroadband.phiPol = normalizedBroadbandMapPolAngles[0];
		normalizedBroadband.thetaPol = -normalizedBroadbandMapPolAngles[1];
		delete normalizedBroadbandMapPol;
//		delete gROOT -> FindObject("peakNormalizedResponseMap");
//		TH2D normalizedBroadbandMapUnpol = makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol, true);
//		normalizedBroadband.peakUnpol = normalizedBroadbandMapUnpol.GetMaximum();
//		vector<double> normalizedBroadbandMapUnpolAngles = getInterferometricPeakLocation(& normalizedBroadbandMapUnpol);
		TH2D * normalizedBroadbandMapUnpol = new TH2D(makePeakInterferometricMap(eventNumber, phiRough, -thetaRough, AnitaPol::kNotAPol, true));
		normalizedBroadband.peakUnpol = normalizedBroadbandMapUnpol -> GetMaximum();
		vector<double> normalizedBroadbandMapUnpolAngles = getInterferometricPeakLocation(normalizedBroadbandMapUnpol);
		normalizedBroadband.phiUnpol = normalizedBroadbandMapUnpolAngles[0];
		normalizedBroadband.thetaUnpol = -normalizedBroadbandMapUnpolAngles[1];
		delete normalizedBroadbandMapUnpol;
//		delete gROOT -> FindObject("peakNormalizedResponseMap");

		pointingTree.Fill();

//		//  Let's see if this helps with clean up? (https://root-forum.cern.ch/t/warning-in-troot-append-replacing-existing-th1-radiance-31-potential-memory-leak/20716)
//		gROOT -> Reset();
//
//		TSeqCollection * canvases = gROOT -> GetListOfCanvases();
//		TIter next(gROOT -> GetListOfCanvases());
//		while (c = (TCanvas*) next()) delete c;
	}

	//  Write the TTree to file.
	pointingFile.cd();
	pointingTree.Write();

	//  Close files.
	sumFile.Close();
	pointingFile.Close();
}



int main(int argc, char * argv[]) {

//	ROOT::EnableThreadSafety();

	int part;  //  Input argument.

	if (argc != 2) {

		cerr << "Usage : " << argv[0] << " [part]" << endl;

		return 1;

	} else part = atoi(argv[1]);

	//  Expedite processing and use sine subtract cache.
	UCorrelator::SineSubtractFilter::setUseCache(true);

	//  Recreate file in case it has to be rewritten.
	TFile pointingFile(TString::Format("testMapPointingFiles/testMapPointing-part_%d.root", part), "recreate");
	pointingFile.Close();

	addPolTree(part);
	addPolTree(part, true);

	return 1;
}
