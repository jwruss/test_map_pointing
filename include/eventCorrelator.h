/*
 * eventCorrelator.h
 *
 *  Created on: Jun 29, 2020
 *      Author: John Russell
 *
 *  Header file from which to reference correlator functions.
 */

using namespace std;

#include <cmath>
#include <omp.h>
#include "AnitaDataset.h"
#include "AnitaGeomTool.h"
#include "AntennaPositions.h"
#include "BasicFilters.h"
#include "UCFilters.h"  //  To be able to use the BH13Filter without declaring a "cfg" object.
#include "FilterStrategy.h"
#include "FilteredAnitaEvent.h"
#include "BH13Filter.h"
#include "graphTools.h"
#include "TString.h"
#include "TH2.h"
#include "TStyle.h"
#include "AnitaVersion.h"
#include "TFile.h"
#include "TSpline.h"
#include "TCollection.h"
#include "TKey.h"
//#include "TROOT.h"


/*
 * Return coarsely-binned histogram of time delays between antenna baselines.
 */
TH2F * getDeltaTCoarse(const char * mapName);


/*
 * Return finely-binned histogram of time delays between antenna baselines.
 */
TH2F * getDeltaTFine(const char * mapName);


/*
 * Return coarsely-binned histogram of spherical cosine products between antenna baselines.
 */
TH2F * getSphCosProductCoarse(const char * mapName);


/*
 * Return finely-binned histogram of spherical cosine products between antenna baselines.
 */
TH2F * getSphCosProductFine(const char * mapName);


/*
 * Get antenna pairs used to construct interferometric maps.
 */
vector<pair<int, int>> getAntennaPairs(vector<int> neighboringAntennas = {});


/*
 * Given an event, construct the cross-correlation between two antennas for a given polarization.
 */
void fillPairMap(FilteredAnitaEvent * filtEvent, TH2D * responseMap, pair<int, int> antPair, AnitaPol::AnitaPol_t pol);


/*
 * Function which iteratively calls "fillPairMap()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllMaps(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas = {});


/*
 * For a given event, return a vector containing sums of square voltages for each antenna. Vector has dimensions of
 * [NUM_SEAVEYS][2]. where second index corresponds to 0 for Hpol, 1 for Vpol.
 */
vector<pair<double, double>> getEventTotalPowers(FilteredAnitaEvent * filtEvent);


/*
 * For normalization purposes. For both H and Vpol channel of each antenna, calculate the "total power" histogram.
 * For a pair of antennas at a given azimuth and elevation, normalization is given by taking the product of the two "total power" histograms
 * at that point and taking the square root of that product.
 */
vector<TH2D> getTotalPowerMaps(FilteredAnitaEvent * filtEvent, bool isFine = false);


/*
 * Each antenna cross-correlation has an upper bound given by the Cauchy-Schwarz inequality. The upper bound here is over both horizontal and vertical polarization.
 */
void fillPairCoverageMap(TH2D * coverageMap, vector<TH2D> totalPowerMaps, pair<int, int> antPair);


/*
 * Function which iteratively calls "fillPairCoverageMap()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllCoverageMaps(TH2D * coverageMap, vector<TH2D> totalPowerMaps, vector<int> neighboringAntennas = {});


/*
 * For a given event, calls to "fillAllMaps()" to construct corresponding unnormalized interferoemtric maps.
 */
vector<TH2D> makeUnnormalizedEventInterferometricMaps(int eventNum, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * For a given event, calls to "fillAllMaps()" and "fillAllCoverageMaps()" to construct corresponding interferometric maps.
 */
vector<TH2D> makeEventInterferometricMaps(int eventNum, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * For a given interferometric map, determine the locations in phi and theta of the peak in the map.
 */
vector<double> getInterferometricPeakLocation(TH2D * responseMap);


/*
 * Given where an interferometric peak occurs, determine which phi sector it is, then return antenna numbers up to three phi sectors apart.
 */
vector<int> getNeighboringAntennas(double peakPhi);


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 * No normalization applied.
 */
TH2D makePeakUnnormalizedInterferometricMap(TH2D * responseMap, int eventNum, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakUnnormalizedInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 */
TH2D makePeakInterferometricMap(TH2D * responseMap, int eventNum, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, bool useBroadband = false, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * For purposes of testing noncoplanar baselines, create cross-correlation
 * interferometric maps without curvature correction.
 */
void fillPairFlatMap(FilteredAnitaEvent * filtEvent, TH2D * responseMap, pair<int, int> antPair, AnitaPol::AnitaPol_t pol);


/*
 * Function which iteratively calls "fillMapsFlatPair()" for all relevant antenna pairs, up to two phi sectors apart.
 */
void fillAllFlatMaps(FilteredAnitaEvent * filtEvent, TH2D * responseMap, AnitaPol::AnitaPol_t pol, vector<int> neighboringAntennas = {});


/*
 * Same principal as fillPairCoverageMap(), but no curvature included.
 * But rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillPairFlatCoverageMap(TH2D * coverageMap, vector<vector<double>> totalPowers, pair<int, int> antPair);


/*
 * Function which iteratively calls "fillFlatCoverageMapsPair()" for all relevant antenna pairs, up to two phi sectors apart.
 * Rather than an input of vector<TH2D> totalPowerMaps, use vector<vector<double>> totalPowers.
 */
void fillFlatCoverageMapsAll(TH2D * coverageMap, vector<vector<double>> totalPowers, vector<int> neighboringAntennas = {});


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferoemtric map about the largest peak in the input map.
 * No normalization applied. Unlike makePeakUnnormalizedInterferometricMap(), uses flat maps.
 */
TH2D makePeakUnnormalizedFlatInterferometricMap(TH2D * responseMap, int eventNum, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);

/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakUnnormalizedFlatInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * Using one of the TH2D objects produced in "makeEventInterferometricMaps()", construct a more finely-binned interferometric map about the largest peak in the input map.
 * Unlike makePeakInterferometricMap(), uses flat maps.
 */
TH2D makePeakFlatInterferometricMap(TH2D * responseMap, int eventNum, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);


/*
 * Different version of above, where you have interferometric peak location on hand, but you have to specify polarization of peak interferometric map.
 */
TH2D makePeakFlatInterferometricMap(int eventNum, double peakPhi, double peakNegTheta, AnitaPol::AnitaPol_t pol, TString filterString = "sinsub_10_3_ad_2", int anitaVer = AnitaVersion::get(), int iceMCRun = 0);
