#ifndef TIMEOFFLIGHT_DETECTOR_H
#define TIMEOFFLIGHT_DETECTOR_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>

class detector {
private:
    std::vector<std::vector<int>> densityMatrix;
    double posEmmiters[2048][2];
    std::vector<std::vector<int>> waterToF;
    std::vector<std::vector<int>> experimentalToF;
    bool info;
public:
    detector(const std::vector<std::vector<int>> waterToF,
              const std::vector<std::vector<int>> experimentalToF,
              const bool& info);
    void findingDifferences(const int threshold, const int indent);
    void printInFile(std::string fileName);
    void testPrint(std::string fileName);
};

#endif //TIMEOFFLIGHT_DETECTOR_H