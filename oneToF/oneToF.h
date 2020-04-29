#ifndef ONETOF_H
#define ONETOF_H

#include <cmath>
#include <vector>
#include <limits>

double disInNumberI(std::vector<int> a, const int &window, const int &tickNumber) {
  double expectedValue = 0, dispersion = 0;
  for (int i = tickNumber; i < tickNumber + window; i++)
    expectedValue += a[i];
  expectedValue /= window;
  for (int i = tickNumber; i < tickNumber + window; i++)
    dispersion += (a[i] - expectedValue) * (a[i] - expectedValue);
  return dispersion / window;
}

double nextDis(std::vector<int> a, const int &window, const double &prevDis,
               const int &nextTickNumber) {
    double expectedValue = 0, sigma = (a[nextTickNumber + window - 1] - a[nextTickNumber - 1]) / window;
    for (int i = 0; i < window; i++)
        expectedValue += a[nextTickNumber - 1 + i];
    return prevDis + sigma * ((window - 1) * (a[nextTickNumber + window - 1]) + (window + 1) * a[nextTickNumber - 1]
    -2 * expectedValue) / window;
}

double averageVariance(std::vector<int> a, const int &window, const double &prevDis,
                       const int &tickNumber) {
  double sum = prevDis, next = prevDis;
  for (int i = tickNumber + 1; i < window; i++) {
    next = nextDis(a, window, prevDis, i);
    sum += next;
  }
  return sum / window;
}

int akaike(std::vector<int> a, const int &window, const int &tickNumber) {
  int signalTick = 0;
  double minAkaike = std::numeric_limits<double>::max();
  for (int i = tickNumber - window; i < tickNumber + window; i++) {
    double leftDis, rightDis;
    leftDis = disInNumberI(a, i - tickNumber + window, tickNumber - window);
    rightDis = disInNumberI(a, window + tickNumber - i, i);
    double tempAkaike = (i - tickNumber + window) * std::log(leftDis) +
                        (window + tickNumber - i) * std::log(rightDis);
    if (tempAkaike < minAkaike) {
        minAkaike = tempAkaike;
        signalTick = i;
    }
  }
  return signalTick;
}

std::vector<int> oneToF(std::vector<int> a, const int &numberOfWaves, const int &windowAkaike,
            const double &cutoffLevel) {
  int m = 0;
  bool searchTerm = false;
  std::vector<int> result;
  double TH = cutoffLevel;
  double S = disInNumberI(a, windowAkaike, 0), averageS;
  for (int i = (windowAkaike + 1); i < (3750 - windowAkaike); i++) {
    S = nextDis(a, windowAkaike, S, i);
//    if (searchTerm) {
//      double tempAverageS = averageVariance(
//          a, windowAkaike, disInNumberI(a, windowAkaike, i - windowAkaike),
//          i - windowAkaike);
//      if (tempAverageS < averageS)
//        averageS = tempAverageS;
//      else {
//        TH = 1.5 * averageS;
//        searchTerm = false;
//      }
//    }
    if (!searchTerm && S > TH) {
      result.push_back(akaike(a, windowAkaike, i));
//      i = result[m] + windowAkaike;
//      averageS = averageVariance(
//          a, windowAkaike, disInNumberI(a, windowAkaike, i - windowAkaike),
//          i - windowAkaike);
//        searchTerm = true;
        m++;
    }
    if (m == numberOfWaves)
        break;
  }
  if (m == 0)
      result.push_back(0);
  return result;
}

#endif // ONETOF_H
