#include "CubicSmile.h"
#include "BSAnalytics.h"
#include <iostream>
#include <cmath>
#include <map>
#include <numeric>

#include "Solver/Eigen/Core"
#include "Solver/LBFGSB.h"
#include "helper.h"

using namespace LBFGSpp;
typedef double Scalar;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

class CubicSmileObjective
{
private:
  std::vector<TickData> volTickerSnap;
  double forward;
  double timeToExp;
  // double acceptable_error;
  bool initialized_best_parameters;

public:
  double best_error;
  Eigen::VectorXd best_parameters;

  CubicSmileObjective(const std::vector<TickData> &volTickerSnap_, double fwd, double T)
      : volTickerSnap(volTickerSnap_),
        forward(fwd),
        timeToExp(T),
        // acceptable_error(acceptable_error_),
        initialized_best_parameters(false)
  {
    best_error = std::numeric_limits<double>::max();
  }

  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad)
  {
    // std::cout << "-----------inside optimisation ----------" << std::endl;
    // std::cout << "x = \n"
    //           << x.transpose() << std::endl;
    double atmvol = x(0);
    double bf25 = x(1);
    double rr25 = x(2);
    double bf10 = x(3);
    double rr10 = x(4);

    // Construct CubicSmile object
    CubicSmile cs(forward, timeToExp, atmvol, bf25, rr25, bf10, rr10);
    double error = CalculateFittingError(volTickerSnap, cs);

    // Initialize the best_parameters with the first x if not already initialized
    if (!initialized_best_parameters)
    {
      best_parameters = x;
      best_error = error;
      initialized_best_parameters = true;
    }

    // Compute gradient of objective function with respect to parameters
    const double delta = 0.001; // or an appropriate small number for your problem

    for (int i = 0; i < x.size(); ++i)
    {
      Eigen::VectorXd x_plus_delta = x;
      x_plus_delta(i) += delta;

      // Recompute the parameters based on updated vector
      atmvol = x_plus_delta(0);
      bf25 = x_plus_delta(1);
      rr25 = x_plus_delta(2);
      bf10 = x_plus_delta(3);
      rr10 = x_plus_delta(4);

      // Construct another CubicSmile object
      CubicSmile cs_plus_delta(forward, timeToExp, atmvol, bf25, rr25, bf10, rr10);

      // Compute the error for this altered parameter set
      double error_plus_delta = CalculateFittingError(volTickerSnap, cs_plus_delta);

      // The gradient is the change in error divided by the change in the parameter
      grad(i) = (error_plus_delta - error) / delta;
    }

    // std::cout << "Gradient values: "
    //           << "grad[0] = " << grad[0] << ", "
    //           << "grad[1] = " << grad[1] << ", "
    //           << "grad[2] = " << grad[2] << ", "
    //           << "grad[3] = " << grad[3] << ", "
    //           << "grad[4] = " << grad[4] << ","
    //           << "error = " << error << std::endl;
    // std::cout << "x = \n"
    //           << x.transpose() << std::endl;

    // Check if the error is less than the acceptable error threshold and update the best error and parameter values

    if (error < best_error)
    {
      best_error = error;
      best_parameters = x;
    }

    return error;
  }
};

CubicSmile CubicSmile::FitSmile(const std::vector<TickData> &volTickerSnap)
{
  double fwd, T, atmvol, bf25, rr25, bf10, rr10;
  // TODO (step 3): fit a CubicSmile that is close to the raw tickers
  // - make sure all tickData are on the same expiry and same underlying
  // - get latest underlying price from all tickers based on LastUpdateTimeStamp
  // - get time to expiry T
  // - fit the 5 parameters of the smile, atmvol, bf25, rr25, bf10, and rr10 using L-BFGS-B solver, to the ticker data
  // ....
  // after the fitting, we can return the resulting smile
  // std::cout << volTickerSnap.size() << std::endl;

  std::map<datetime_t, int> expDate;
  for (const TickData &tickData : volTickerSnap)
  {
    expDate[tickData.ExpiryDate]++;
  }

  // for (const auto &entry : expDate)
  // {
  //   std::cout << "Unique date:  " << entry.first << " , count: " << entry.second << std::endl;
  // }

  double undPriceSum = std::accumulate(volTickerSnap.begin(), volTickerSnap.end(), 0.0,
                                       [](double acc, const TickData &data)
                                       {
                                         return acc + data.UnderlyingPrice;
                                       });

  // use average since there are some slight difference of the timestamp which resulted in different underlying price
  fwd = undPriceSum / static_cast<double>(volTickerSnap.size());

  datetime_t currentTimeStamp(volTickerSnap[0].LastUpdateTimeStamp / 1000);
  datetime_t expiryDate = volTickerSnap[0].ExpiryDate;
  // getting the time to maturity
  T = (expiryDate - currentTimeStamp) / 365;
  T = T < 0.001 ? 0.001 : T;

  atmvol = GetATMVolatility(volTickerSnap, fwd);

  // we use quick delta: qd = N(log(F/K / (atmvol) / sqrt(T))
  double stdev = atmvol * sqrt(T);
  double k_qd90 = quickDeltaToStrike(0.9, fwd, stdev);
  double k_qd75 = quickDeltaToStrike(0.75, fwd, stdev);
  double k_qd25 = quickDeltaToStrike(0.25, fwd, stdev);
  double k_qd10 = quickDeltaToStrike(0.10, fwd, stdev);

  double vol90 = interpolateQuickDeltaIV(volTickerSnap, k_qd90, false);
  double vol75 = interpolateQuickDeltaIV(volTickerSnap, k_qd75, false);
  double vol25 = interpolateQuickDeltaIV(volTickerSnap, k_qd25, true);
  double vol10 = interpolateQuickDeltaIV(volTickerSnap, k_qd10, true);

  bf25 = ((vol75 + vol25) / 2.0) - atmvol;
  bf10 = ((vol90 + vol10) / 2.0) - atmvol;
  rr25 = vol25 - vol75;
  rr10 = vol10 - vol90;

  auto sm = CubicSmile(fwd, T, atmvol, bf25, rr25, bf10, rr10);
  std::tie(sm.maxOpenInterest, sm.maxSpread) = GetMaxValues(volTickerSnap);
  double initialError = CalculateFittingError(volTickerSnap, sm);

  // Define parameters
  // const int n = 5;
  LBFGSBParam<Scalar> param;
  // min_step and max_step set the lower and upper bounds on the step size that the line search algorithm can take in any single iteration. If you want to effectively lower the learning rate, you could decrease the max_step value.
  // param.max_step = 0.01; // Lower max_step value

  // ftol and wolfe control the accuracy of the line search routine. Lowering ftol or increasing wolfe can make the line search more stringent, which might effectively result in smaller steps, i.e., lower learning rate.
  // param.ftol = 1e-6;  // Lower ftol value
  // param.wolfe = 0.95; // Increase wolfe value
  param.max_linesearch = 100;

  // Declare and initialize the acceptable error variable
  double acceptable_error = 0.0001;
  // // Set the acceptable error in the LBFGSBParam object
  // param.epsilon = acceptable_error;

  LBFGSBSolver<Scalar> solver(param);
  // OptimiserFunctionObj fun(volTickerSnap, fwd, T);

  // Initial guess for parameters
  Eigen::VectorXd x0(5);
  // Variable bounds
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(5, -10);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(5, 10);
  // Initial values
  // Vector x = Vector::Constant(n, 0.0);
  x0 << atmvol, bf25, rr25, bf10, rr10;

  // std::cout << "-----------Before optimisation ----------" << std::endl;
  // std::cout << "x = \n"
  //           << x0.transpose() << std::endl;

  CubicSmileObjective fun(volTickerSnap, fwd, T);

  Scalar fx;
  try
  {
    int niter = solver.minimize(fun, x0, fx, lb, ub);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n"
              << x0.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    // Check if the final error exceeds the acceptable threshold
    if (fx > acceptable_error)
    {
      std::cout << "Optimization stopped due to error exceeding the acceptable threshold." << std::endl;
      // Optionally, you can access the best error and parameter values seen so far using the `fun` object.
      std::cout << "Best error: " << fun.best_error << std::endl;
      std::cout << "Best parameters: " << fun.best_parameters.transpose() << std::endl;
    }
    else
    {
      std::cout << "Optimization completed successfully." << std::endl;
      std::cout << "Final error: " << fx << std::endl;
      std::cout << "Final parameters: " << x0.transpose() << std::endl;
    }
  }
  catch (const std::runtime_error &e)
  {
    std::cout << "-----------After optimisation run time limit----------" << std::endl;
    std::cout << "x = \n"
              << x0.transpose() << std::endl;
    std::cerr << "Runtime error: " << e.what() << std::endl;

    x0 = fun.best_parameters;
    std::cout << "Using best parameters" << std::endl;
    std::cout << "x = \n"
              << x0.transpose() << std::endl;
    return CubicSmile(fwd, T, x0(0), x0(1), x0(2), x0(3), x0(4), {atmvol, bf25, rr25, bf10, rr10}, initialError);
  }
  catch (...)
  {
    std::cerr << "-----------An unknown error occurred.-----------" << std::endl;
    std::cout << "x = \n"
              << x0.transpose() << std::endl;
    x0 = fun.best_parameters;
    std::cout << "Using best parameters" << std::endl;
    std::cout << "x = \n"
              << x0.transpose() << std::endl;
    return CubicSmile(fwd, T, x0(0), x0(1), x0(2), x0(3), x0(4), {atmvol, bf25, rr25, bf10, rr10}, initialError);
  }

  return CubicSmile(fwd, T, x0(0), x0(1), x0(2), x0(3), x0(4), {atmvol, bf25, rr25, bf10, rr10}, initialError);
}

CubicSmile::CubicSmile(double underlyingPrice, double T, double atmvol,
                       double bf25, double rr25, double bf10, double rr10, std::vector<double> init_guess, double init_error)
{
  // std::cout << "constructing Cubic Smile" << std::endl;
  // convert delta marks to strike vol marks, setup strikeMarks, then call BUildInterp
  double v_qd90 = atmvol + bf10 - rr10 / 2.0;
  double v_qd75 = atmvol + bf25 - rr25 / 2.0;
  double v_qd25 = atmvol + bf25 + rr25 / 2.0;
  double v_qd10 = atmvol + bf10 + rr10 / 2.0;

  // we use quick delta: qd = N(log(F/K / (atmvol) / sqrt(T))
  double stdev = atmvol * sqrt(T);
  double k_qd90 = quickDeltaToStrike(0.9, underlyingPrice, stdev);
  double k_qd75 = quickDeltaToStrike(0.75, underlyingPrice, stdev);
  double k_qd25 = quickDeltaToStrike(0.25, underlyingPrice, stdev);
  double k_qd10 = quickDeltaToStrike(0.10, underlyingPrice, stdev);

  strikeMarks.push_back(std::pair<double, double>(k_qd90, v_qd90));
  strikeMarks.push_back(std::pair<double, double>(k_qd75, v_qd75));
  strikeMarks.push_back(std::pair<double, double>(underlyingPrice, atmvol));
  strikeMarks.push_back(std::pair<double, double>(k_qd25, v_qd25));
  strikeMarks.push_back(std::pair<double, double>(k_qd10, v_qd10));

  future_price = underlyingPrice;
  primer_guess = init_guess;
  primer_error = init_error;

  BuildInterp();
  // BuildInterpNotAKnot();
}

vector<pair<double, double>> CubicSmile::GetStrikeMarks() const
{
  return strikeMarks;
}

void CubicSmile::BuildInterp()
{
  int n = strikeMarks.size();
  // end y' are zero, flat extrapolation
  double yp1 = 0;
  double ypn = 0;
  y2.resize(n);
  vector<double> u(n - 1);

  y2[0] = -0.5;
  u[0] = (3.0 / (strikeMarks[1].first - strikeMarks[0].first)) *
         ((strikeMarks[1].second - strikeMarks[0].second) / (strikeMarks[1].first - strikeMarks[0].first) - yp1);

  // Calculate second derivatives at each point
  for (int i = 1; i < n - 1; i++)
  {
    double sig = (strikeMarks[i].first - strikeMarks[i - 1].first) / (strikeMarks[i + 1].first - strikeMarks[i - 1].first);
    double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (strikeMarks[i + 1].second - strikeMarks[i].second) / (strikeMarks[i + 1].first - strikeMarks[i].first) - (strikeMarks[i].second - strikeMarks[i - 1].second) / (strikeMarks[i].first - strikeMarks[i - 1].first);
    u[i] = (6.0 * u[i] / (strikeMarks[i + 1].first - strikeMarks[i - 1].first) - sig * u[i - 1]) / p;
  }

  double qn = 0.5;
  double un = (3.0 / (strikeMarks[n - 1].first - strikeMarks[n - 2].first)) *
              (ypn - (strikeMarks[n - 1].second - strikeMarks[n - 2].second) / (strikeMarks[n - 1].first - strikeMarks[n - 2].first));

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  //  std::cout << "y2[" << n-1 << "] = " << y2[n-1] << std::endl;
  for (int i = n - 2; i >= 0; i--)
  {
    y2[i] = y2[i] * y2[i + 1] + u[i];
    //    std::cout << "y2[" << i << "] = " << y2[i] << std::endl;
  }
}

/**
 * In the revised Vol() function, when we need to extrapolate to the left (when strike is less than the smallest strike in strikeMarks),
 * we pretend that the strike lies within the first interval of our strikeMarks (between strikeMarks[0] and strikeMarks[1]).
 * We then use the same cubic spline that was built for that interval to estimate the volatility for the strike. This is done by setting i = 1 if i == 0.
 * Similarly, when we need to extrapolate to the right (when strike is greater than the largest strike in strikeMarks),
 * we pretend that the strike lies within the last interval of our strikeMarks (between strikeMarks[n-2] and strikeMarks[n-1],
 * where n is the size of strikeMarks). We then use the same cubic spline that was built for that interval to estimate the volatility for the strike.
 * This is done by setting i = strikeMarks.size() - 1 if i == strikeMarks.size().
 * Once we have the bracketing interval for the strike, whether it's due to actual interpolation or pretended for extrapolation, we use the cubic spline formula to estimate the volatility.
 * The cubic spline formula consists of a weighted average of the volatilities at the two ends of the bracketing interval (strikeMarks[i-1].second and strikeMarks[i].second),
 * and two terms involving the second derivatives at the two ends (y2[i-1] and y2[i]) which ensure that the function is smooth and has continuous first and second derivatives.
 * While this method can provide reasonable results when the strike is not too far outside the range of strikeMarks, it is still a form of extrapolation,
 * and its accuracy will decrease the further strike is from the range of strikeMarks.
 * The extrapolation is essentially assuming that the trend observed within the range (as characterized by the cubic spline) continues outside the range.
 * If this assumption is not accurate, then the extrapolated volatilities may not be accurate.
 */
double CubicSmile::Vol(double strike) const
{
  unsigned i;
  // we use trivial search, but can consider binary search for better performance
  for (i = 0; i < strikeMarks.size(); i++)
    if (strike <= strikeMarks[i].first)
      break; // i stores the index of the right end of the bracket

  // // interpolate or extrapolate
  // if (i == 0) // extrapolation to the left
  //   i = 1;
  // if (i == strikeMarks.size()) // extrapolation to the right
  //   i = strikeMarks.size() - 1;

  // extrapolation
  // this one return vol at the smallest and largest strike no extrapolation
  if (i == 0)
    return strikeMarks[i].second;
  if (i == strikeMarks.size())
    return strikeMarks[i - 1].second;

  // Now we are sure that strike is bracketed by strikeMarks[i-1] and strikeMarks[i]

  // interpolate or extrapolate using the cubic spline
  double h = strikeMarks[i].first - strikeMarks[i - 1].first;
  double a = (strikeMarks[i].first - strike) / h;
  double b = 1 - a;
  double c = (a * a * a - a) * h * h / 6.0;
  double d = (b * b * b - b) * h * h / 6.0;
  return a * strikeMarks[i - 1].second + b * strikeMarks[i].second + c * y2[i - 1] + d * y2[i];
}

/**
 * The fitting error some how become very bad shall not use this
 */
void CubicSmile::BuildInterpNotAKnot()
{
  int n = strikeMarks.size();
  y2.resize(n);
  vector<double> u(n - 1);

  double sig, p;

  y2[0] = 0.0;
  u[0] = (3.0 / (strikeMarks[1].first - strikeMarks[0].first)) *
         ((strikeMarks[1].second - strikeMarks[0].second) / (strikeMarks[1].first - strikeMarks[0].first));

  for (int i = 1; i < n - 1; i++)
  {
    sig = (strikeMarks[i].first - strikeMarks[i - 1].first) / (strikeMarks[i + 1].first - strikeMarks[i - 1].first);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (strikeMarks[i + 1].second - strikeMarks[i].second) / (strikeMarks[i + 1].first - strikeMarks[i].first) - (strikeMarks[i].second - strikeMarks[i - 1].second) / (strikeMarks[i].first - strikeMarks[i - 1].first);
    u[i] = (6.0 * u[i] / (strikeMarks[i + 1].first - strikeMarks[i - 1].first) - sig * u[i - 1]) / p;
  }

  y2[n - 1] = 0.0;

  y2[n - 1] = (y2[n - 3] / 4.0 - y2[n - 2] / 2.0 - u[n - 2]) / (0.5 - y2[n - 2]);
  for (int i = n - 2; i >= 0; i--)
  {
    y2[i] = y2[i] * y2[i + 1] + u[i];
  }

  y2[0] = y2[2] / 4.0 - y2[1] / 2.0;
}