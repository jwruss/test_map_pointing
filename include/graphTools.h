using namespace std;

#include <cmath>
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "FFTWComplex.h"
#include "FFTtools.h"
#include "TAxis.h"

/*  TGraph::Integral() returns the area within a closed polygon of a TGraph after bin reshuffling,
 *  which is a different operation than integrating the function corresponding to TGraph::GetY() over
 *  the domain corresponding to TGraph::GetX(). Using trapezoidal rule and accounting for potential uneven
 *  sampling in the TGraph, we explicitly calculate the integral here.
 */
double getGraphIntegral(const TGraph * grInPtr);


/*  The weighting done in FFTtools::getCorrelation() appears to be off.
 */
TGraph getCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr);


/*  To use interpolation of TGraphs used in getCorrGraph() instead.
 */
TGraph getInterpolatedCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr, double deltaT = 0.1);


/*  Normalized version of getCorrGraph().
 *  We correct it here. Dropping "normalization" from the name of getCorrGraph() comes from the fact
 *  that what FFTtools::getCorrelationGraph() actually attempts to produce is the cross-covariance
 *  graph. Its normalization is what is usually referred to as the correlation, or cross-correlation.
 */
TGraph getNormalizedCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr);


/*  To use interpolation of TGraphs used in getNormalizedCorrGraph() instead.
 */
TGraph getInterpolatedNormalizedCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr, double deltaT = 0.1);


/*  Create a truncation of input TGraph grInPtr to a ns-wide window centered about the pulse's zero crossing between maximum and minimum values.
 *  The input winFullWidth is the ns-width over the full window; it's symmetrically reduced when outside the scope of grInPtr's range.
 */
TGraph getWindowedGraph(const TGraph * grInPtr, double winFullWidth = 5);


/*  To use interpolation of TGraph used in getWindowedGraph() instead.
 */
TGraph getWindowedInterpolatedGraph(const TGraph * grInPtr, double deltaT = 0.1, double winFullWidth = 5);


/*  Get TGraph which is an affine transformation of input TGraph.
 *  Setting "scale = -1, offset = 0" flips polarity.
 */
TGraph getAffineTransformedGraph(const TGraph * grInPtr, double scale = 1, double offset = 0);


/*  Take an input TGraph and generate another TGraph corresponding to its numeric central derivative.
 *  Technically, the input TGraph doesn't have to be uniformly sampled.
 */
TGraph getDerivGraph(const TGraph * grInPtr);


/*  Same as getDerivGraph(), but output TGraph's y-component has polarity reversed.
 */
TGraph getNegDerivGraph(const TGraph * grInPtr);


/*  Take an input TGraph and generate another TGraph corresponding to its numerical antiderivative.
 *  Technically, the input TGraph doesn't have to be uniformly sampled. Taking grInPtr -> GetY()
 *  and integrating up from the lower bound element, down from the upper bound element, and then
 *  averaging these results to find an accurate and invertible antiderivative, calling the method
 *  employed here a "superposition method" seems appropriate.
 */
TGraph getAntiderivGraph(const TGraph * grInPtr);


/*  In principal this computes the same as "getAntiDerivGraph()", but the method here relies
 *  on midpoint integration method. Should be faster to compute as well, and work best for
 *  uniformly sampled data.
 */
TGraph getMidpointAntiderivGraph(const TGraph * grInPtr);


/*  In principal this computes the same as "getAntiDerivGraph()" and "getAntiDerivMidpointGraph()",
 *  but the method here relies on trapezoidal integration method. Should be faster to compute as well.
 */
TGraph getTrapezoidalAntiderivGraph(const TGraph * grInPtr);


/*  Modifying scipy.fftpack.diff from SciPy
 *  <https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.fftpack.diff.html>
 *  for usage with FFTtools. Assumes the waveform is uniformly sampled, zero-meaned when
 *  order is negative, isn't constant after differentiation, and the result is expected to be
 *  real over the input domain, so don't use it otherwise.
 */
TGraph getIFFTDiffGraph(const TGraph * grInPtr, double order = 1, int branchOrder = 0);


/*  Same as getNegDerivGraph(), except using getIFFTDiffGraph() instead.
 */
TGraph getIFFTNegDerivGraph(const TGraph * grInPtr);
