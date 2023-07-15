#ifndef HELPER_H
#define HELPER_H

#include <utility>
#include <vector>
#include <chrono>
#include "Msg.h"
#include "CubicSmile.h"

#include "helper.h"

inline std::pair<double, double> GetMaxValues(const std::vector<TickData> &volTickerSnap)
{
  double maxOpenInterest = 0.0;
  double maxSpread = 0.0;

  for (const TickData &tickData : volTickerSnap)
  {
    // Update the maximum open interest if necessary
    if (tickData.OpenInterest > maxOpenInterest)
      maxOpenInterest = tickData.OpenInterest;

    // Calculate the spread
    double spread = tickData.BestAskPrice - tickData.BestBidPrice;

    // Update the maximum spread if necessary
    if (spread > maxSpread)
      maxSpread = spread;
  }

  return std::make_pair(maxOpenInterest, maxSpread);
}

inline double CalculateWeight(const TickData &tickdata, double maxOpenInterest, double maxSpread)
{

  double normOpenInterest = tickdata.OpenInterest / maxOpenInterest;
  double normSpread = (tickdata.BestAskPrice - tickdata.BestBidPrice) / maxSpread;
  // std::cout << "normOpenInterest = " << normOpenInterest << std::endl;
  // std::cout << "normSpread = " << normSpread << std::endl;
  return (normOpenInterest + normSpread) / 2.0;
}

inline double CalculateFittingError(const std::vector<TickData> &volTickerSnap, const CubicSmile &sm)
{
  double fittingError = 0.0;
  double sumWeights = 0.0;

  // GetMaxValues(volTickerSnap); //

  // Get the maximum values for weight calculation
  std::pair<double, double> maxValues = GetMaxValues(volTickerSnap);
  double maxOpenInterest = maxValues.first;
  double maxSpread = maxValues.second;

  for (const TickData &tickData : volTickerSnap)
  {
    // Calculate weight based on liquidity, open interest, bid-ask spread
    double weight = CalculateWeight(tickData, maxOpenInterest, maxSpread);
    // std::cout << "TickData contract = " << tickData.ContractName << std::endl;
    // std::cout << "weight = " << weight << std::endl;
    // Calculate average implied volatility
    double sigma_i = (tickData.BestBidIV + tickData.BestAskIV) / 2.0;

    // Calculate model implied volatility
    // TODO: Check why is this nan
    double sigma_ki = sm.Vol(tickData.Strike);

    // Calculate difference and multiply by weight
    double diff = std::abs(sigma_i - sigma_ki);
    double weightedDiff = diff * weight;
    // std::cout << "sigma difference = " << diff << std::endl;
    // std::cout << "weightedDiff = " << weightedDiff << std::endl;

    // Update fitting error and sum of weights
    fittingError += weightedDiff;
    sumWeights += weight;
  }

  // Divide the sum of weighted differences by the sum of weights to get the fitting error
  if (sumWeights != 0.0)
  {
    fittingError /= sumWeights;
  }

  // std::cout << "fittingError = " << fittingError << std::endl;
  // std::cout << "sumWeights = " << sumWeights << std::endl;
  return fittingError;
}

// function to convert epoch time in millisec to string eg 2022-07-02T01:38:07.232Z
inline std::string convert_msec_to_utc_string(uint64_t msec)
{

  // std::cout << "incoming msec = " << msec << std::endl;

  datetime_t datetime = datetime_t(msec / 1000); // Convert milliseconds to seconds and create datetime_t object

  // Create a formatted string in UTC format
  char buffer[30];
  std::sprintf(buffer, "%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",
               datetime.year, datetime.month, datetime.day,
               datetime.hour, datetime.min, datetime.sec, (int)(msec % 1000));
  // std::cout << "output string = " << std::string(buffer) << std::endl;
  return std::string(buffer);
}

inline std::string getMonthAbbreviation(int month)
{
  switch (month)
  {
  case 1:
    return "JAN";
  case 2:
    return "FEB";
  case 3:
    return "MAR";
  case 4:
    return "APR";
  case 5:
    return "MAY";
  case 6:
    return "JUN";
  case 7:
    return "JUL";
  case 8:
    return "AUG";
  case 9:
    return "SEP";
  case 10:
    return "OCT";
  case 11:
    return "NOV";
  case 12:
    return "DEC";
  default:
    return "XXX";
  }
}

#endif