#ifndef _CUBICSMILE_H
#define _CUBICSMILE_H

#include <vector>
#include <utility>
#include "Msg.h"

using namespace std;

// CubicSpline interpolated smile, extrapolate flat
class CubicSmile
{
public:
  // FitSmile creates a Smile by fitting the smile params to the input tick data, it assume the tickData are of the same expiry
  static CubicSmile FitSmile(const std::vector<TickData> &);
  // constructor, given the underlying price and marks, convert them to strike to vol pairs (strikeMarks), and construct cubic smile
  CubicSmile(double underlyingPrice, double T, double atmvol, double bf25, double rr25, double bf10, double rr10, vector<double> init_guess = {0, 0, 0, 0, 0}, double init_error = 0.0); // convert parameters to strikeMarks, then call BuildInterp() to create the cubic spline interpolator
  double Vol(double strike) const;

  vector<pair<double, double>> GetStrikeMarks() const;
  double future_price;
  vector<double> primer_guess;
  double primer_error; // interpolate
  double maxOpenInterest, maxSpread;

private:
  void BuildInterp();
  void BuildInterpNotAKnot();
  // strike to implied vol marks
  vector<pair<double, double>> strikeMarks;
  vector<double> y2; // second derivatives

};

#endif
