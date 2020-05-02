#include "detector.h"
#include <iostream>


detector::detector(const std::vector<std::vector<int>> waterToF,
                   const std::vector<std::vector<int>> experimentalToF,
                   const bool& info)
                   : waterToF(waterToF), experimentalToF(experimentalToF),
                   info(info) {
    int j = 0;
    double r = 1024;
    for (int i = 0; i < 2048 + 56; i ++) {
        if (!((i % (256 + 8)) < 8  && i > 8)) {
            posEmmiters[j][0] = r + r * std::cos(2*M_PI * i / (2048 + 64));
            posEmmiters[j][1] = r + r * std::sin(2*M_PI * i / (2048 + 64));
            j++;
        }
    }
    densityMatrix = std::vector<std::vector<int>>(2048);
    for (int i = 0; i < 2048; i++)
        densityMatrix[i] = std::vector<int>(2048);
}

void detector::findingDifferences(const int threshold, const int indent) {
    double eps = 0.5;
    int test = 0;
    std::set<int> brokenSensor = {528, 1024, 1025, 1026, 1027, 1028, 1029, 1031, 1152, 1154, 1155, 1156, 1157, 1158, 1188};
    for (int i = 0; i < 2048; i ++) {
        int count = 0;
        for (int j = (i + 255) % 2048; count < 2048 - 2 * indent; j = (j + 1) % 2048) {
            count++;
            if (brokenSensor.count(i) || brokenSensor.count(j))
                continue;
            if (std::abs(waterToF[i][j] - experimentalToF[i][j]) < threshold) {
                test++;
                //line x = x[i]
                if (std::abs(posEmmiters[i][0] - posEmmiters[j][0]) < eps) {
                    for (int x = 1; x < 2047; x++) {
                        if (std::abs(x - (posEmmiters[i][0] + posEmmiters[j][0]) / 2) < eps) {
                            for (int y = 1; y < 2047; y++) {
                                densityMatrix[x][y] ++;
                            }
                            continue;
                        }
                    }
                    continue;
                }
                //line y = y[i]
                else if (std::abs(posEmmiters[i][1] - posEmmiters[j][1]) < eps) {
                    for (int y = 0; y < 2048; y++) {
                        if (std::abs(y - (posEmmiters[i][1] + posEmmiters[j][1]) / 2) < eps) {
                            for (int x = 0; x < 2048; x++) {
                                densityMatrix[x][y] ++;
                            }
                            continue;
                        }
                    }
                    continue;
                }
                //line y = ax + b
                else {
                    bool tempMatrix[2048][2048];
                    for (int x = 0; x < 2048; x++) {
                        for (int y = 0; y < 2048; y++)
                            tempMatrix[x][y] = false;
                    }
                    double a = ((double)(posEmmiters[j][1] - posEmmiters[i][1])/(posEmmiters[j][0] - posEmmiters[i][0])),
                           b = posEmmiters[i][1] - a * posEmmiters[i][0];
                    if (std::abs(a) > 1) {
                        for (int y = 1; y < 2047; y++) {
                            int tempX = (int) ((y + 0.5 - b) / a + 0.5);
                            if (tempX < 2048 && tempX > 0) {
                                tempMatrix[tempX][y] = true;
                                tempMatrix[tempX][y + 1] = true;
                            }
                        }
                    }
                    if (std::abs(a) < 1) {
                        for (int x = 1; x < 2047; x++) {
                            int tempY = (int) (a * (x + 0.5) + b);
                            if (tempY < 2048 && tempY > 0) {
                                tempMatrix[x][tempY] = true;
                                tempMatrix[x + 1][tempY] = true;
                            }
                        }
                    }
                    for (int x = 0; x < 2048; x++) {
                        for (int y = 0; y < 2048; y++) {
                            if (tempMatrix[x][y]) {
                                densityMatrix[x][y] ++;
                            }
                        }
                    }
                }
            }
        }
        if (info)
            std::cout << "\t emmiter number " << i << " is ready!\n";
    }
}

void detector::printInFile(std::string fileName) {
    std::ofstream outPut(fileName);
    for (int i = 0; i < 2048; i++) {
        for (int j = 0; j < 2048; j++) {
            if ((i - 1024) * (i - 1024) + (j - 1024) * (j - 1024) < 1024 * 1024)
                outPut << densityMatrix[i][j] << " ";
            else
                outPut << 0 << " ";
        }
        outPut << "\n";
    }
}
//                std::cout << j << " ";
void detector::testPrint(std::string fileName) {
    std::ofstream outPut(fileName);
    for (int i = 0; i < 2048; i++)
        outPut << posEmmiters[i][0] << " " << posEmmiters[i][1] << "\n";
}
