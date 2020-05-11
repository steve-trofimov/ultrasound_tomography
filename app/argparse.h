//
// Created by anna on 28.04.2020.
//
#ifndef TIMEOFFLIGHT_ARGPARSE_H
#define TIMEOFFLIGHT_ARGPARSE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <boost/program_options.hpp>

#include "read_data.h"
#include "oneToF.h"
#include "detector.h"
#include "rungekutta.h"
#include "systemEquations.h"

typedef std::chrono::high_resolution_clock timer;

//parametrs
const double waterSpeed = 1462.738;
const int size = 2048;
const int tickNumber = 3750;
const int numberOfWaves = 1;
const int signalDuration = 80;
const int delay = 80;
const int windowAkaike = 60;
const double cutoffLevel = 50;
bool timeTest = false;
const double endT = 3., startT = 0., h = 1./1000.;

namespace po = boost::program_options;

void objSpeedOptimizationFromSmall(const po::variables_map &vm) {
    std::string input_ToF, input_detection;
    bool info;
    if (vm.count("input_ToF"))
        input_ToF = vm["input_ToF"].as<std::string>();
    if (vm.count("input_detection"))
        input_detection = vm["input_detection"].as<std::string>();
    if (vm.count("info"))
        info = vm["info"].as<bool>();

    //init emmiters
    std::vector<emmiter> emmiters;

    for (int i = 0; i < 2048 + 56; i++) { // i < 2048 + 56
        if (!((i % (256 + 8)) < 8  && i > 8)) {
            emmiter temp(size / 2 + size / 2 * std::cos(2 * M_PI * i / (size + 64)),
                         size / 2 + size / 2 * std::sin(2 * M_PI * i / (size + 64)),
                         size / 2);
            emmiters.push_back(temp);
        }
    }

    //init ToF
    std::ifstream ToFfile(input_ToF);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ToFMat;
    ToFMat.resize(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int tempToF;
            ToFfile >> tempToF;
            ToFMat(i, j) = tempToF;
        }
    }

    std::vector<double> optSpeedValue;

    for (int i = emmiters.size() - 128; i != 128; i = (i + 16) % emmiters.size()) {

        double minError = std::numeric_limits<double>::max(), resultObjSpeed;

        for(int objSpeed = waterSpeed - 200; objSpeed <= waterSpeed + 200; objSpeed += 10) {

            std::fstream inputMatrix(input_detection);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> fieldMat;
            fieldMat.resize(size, size);

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    int temp;
                    inputMatrix >> temp;
                    if (temp > 300 || ((i - 1024) * (i - 1024) + (j - 1024) * (j - 1024) > 1000 * 1000))
                        fieldMat(i, j) = waterSpeed;
                    else
                        fieldMat(i, j) = objSpeed;
                }
            }

            RungeKutta RK;
            systemEquation numericalRhs(fieldMat, waterSpeed);

            auto startTimer = timer::now();
            double timeForRK = 0;
            auto path = RK.functor(startT, endT, h, emmiters[i], numericalRhs, timeForRK);
            int number = emitterNumberSearch(emmiters, path.back());

            double error = (timeForRK * 22 / 100 / 2048 - (ToFMat(i, number) + 248) / 25 / 1000000) *
                           (timeForRK * 22 / 100 / 2048 - (ToFMat(i, number) + 248) / 25 / 1000000);

            std::chrono::duration<double> dur = timer::now() - startTimer;

            if (info)
                std::cout << "\t Object speed: " << objSpeed << ", emmiter " << i << " is complite, RK time: " << dur.count() << "\n";

            if (error < minError) {
                minError = error;
                resultObjSpeed = objSpeed;
            }
        }
        optSpeedValue.push_back(resultObjSpeed);
        std::cout << "\t Optimal object speed: " << resultObjSpeed << ", for emmiter " << i << "\n";
    }

    double averageResult = 0;
    for (auto iter : optSpeedValue)
        averageResult += iter;
    std::cout << "\t Average result: " << averageResult << "\n";
}

void objSpeedOptimizationFromBig(const po::variables_map &vm) {
    std::string input_ToF, input_detection;
    bool info;
    if (vm.count("input_ToF"))
        input_ToF = vm["input_ToF"].as<std::string>();
    if (vm.count("input_detection"))
        input_detection = vm["input_detection"].as<std::string>();
    if (vm.count("info"))
        info = vm["info"].as<bool>();

    //init emmiters
    std::vector<emmiter> emmiters;

    for (int i = 0; i < 2048 + 56; i++) { // i < 2048 + 56
        if (!((i % (256 + 8)) < 8  && i > 8)) {
            emmiter temp(size / 2 + size / 2 * std::cos(2 * M_PI * i / (size + 64)),
                         size / 2 + size / 2 * std::sin(2 * M_PI * i / (size + 64)),
                         size / 2);
            emmiters.push_back(temp);
        }
    }

    //init ToF
    std::ifstream ToFfile(input_ToF);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ToFMat;
    ToFMat.resize(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int tempToF;
            ToFfile >> tempToF;
            ToFMat(i, j) = tempToF;
        }
    }

    double minError = std::numeric_limits<double>::max(), resultObjSpeed;

    for(int objSpeed = waterSpeed - 20; objSpeed <= waterSpeed + 20; objSpeed += 5) {

        std::fstream inputMatrix(input_detection);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> fieldMat;
        fieldMat.resize(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int temp;
                inputMatrix >> temp;
                if (temp > 300 || ((i - 1024) * (i - 1024) + (j - 1024) * (j - 1024) > 1000 * 1000))
                    fieldMat(i, j) = waterSpeed;
                else
                    fieldMat(i, j) = objSpeed;
            }
        }

        RungeKutta RK;
        systemEquation numericalRhs(fieldMat, waterSpeed);

        double error = 0;

        for (int i = emmiters.size() - 128; i != 128; i = (i + 1) % emmiters.size()) {
            auto startTimer = timer::now();
            double timeForRK = 0;
            auto path = RK.functor(startT, endT, h, emmiters[i], numericalRhs, timeForRK);
            int number = emitterNumberSearch(emmiters, path.back());
            error += (timeForRK * 22 / 100 / 2048 - (ToFMat(i, number) + 248) / 25 / 1000000) *
                     (timeForRK * 22 / 100 / 2048 - (ToFMat(i, number) + 248) / 25 / 1000000);
            std::chrono::duration<double> dur = timer::now() - startTimer;
            if (info)
                std::cout << "\t Object speed: " << objSpeed << ", emmiter " << i << " is complite, RK time: " << dur.count() << "\n";
        }

        std::cout << "\t Object speed: " << objSpeed << " is complite, error: " << error << "\n";

        if (error < minError) {
            minError = error;
            resultObjSpeed = objSpeed;
        }
    }

    std::cout << "\t Result object speed: " << resultObjSpeed << " , result error: " << minError << "\n";
}

void ToF(const po::variables_map &vm) {
    std::string input_dir, output;
    bool info;
    if (vm.count("input_dir"))
        input_dir = vm["input_dir"].as<std::string>();
    if (vm.count("output"))
        output = vm["output"].as<std::string>();
    if (vm.count("info"))
        info = vm["info"].as<bool>();

    auto startTimer = timer::now();

    std::ofstream outAll(output);
    for (int i = 0; i < size; i++) {
        std::vector<std::vector<int>> oneEmmiterPtr = readEmitter(input_dir, i + 1);
        for (int j = 0; j < size; j++) {
            int resultOne = oneToF(oneEmmiterPtr[j], numberOfWaves, windowAkaike, cutoffLevel)[0];
            outAll << resultOne << " ";
        }
        outAll << "\n";
        if (info)
            std::cout << "\t time of flight calculate: " << i + 1 << std::endl;
    }

    if (info) {
        std::chrono::duration<double> dur = timer::now() - startTimer;
        std::cout << "\t all time : " << dur.count() << std::endl;
    }
}

void oneToF(const po::variables_map &vm) {
    std::string input_dir, outputToF, outputRead;
    int numberEmmiter;
    if (vm.count("input_dir"))
        input_dir = vm["input_dir"].as<std::string>();
    if (vm.count("number_emmiter"))
        numberEmmiter = vm["number_emmiter"].as<int>();

    outputToF = "ToF_result_of_emmiter_" + std::to_string(numberEmmiter) + ".txt";
//    outputRead = "read_data_from_emmiter_" + std::to_string(numberEmmiter) + ".txt";

    std::ofstream out(outputToF); //outRead(outputRead);
    std::vector<std::vector<int>> oneEmmiterPtr = readEmitter(input_dir, numberEmmiter + 1);
    for (int j = 0; j < size; j++) {
        int resultOne = oneToF(oneEmmiterPtr[j], numberOfWaves, windowAkaike, cutoffLevel)[0];
        out << resultOne << " ";
//        for (int i = 0; i < oneEmmiterPtr[j].size(); i++) {
//            outRead << oneEmmiterPtr[j][i] << " ";
//        }
//        outRead << "\n";
    }
    std::cout << "Done!\n";
}

void detection(const po::variables_map &vm) {
    std::string input_water, input_exp, output;
    bool info;
    if (vm.count("input_water"))
        input_water = vm["input_water"].as<std::string>();
    if (vm.count("input_exp"))
        input_exp = vm["input_exp"].as<std::string>();
    if (vm.count("output"))
        output = vm["output"].as<std::string>();
    if (vm.count("info"))
        info = vm["info"].as<bool>();

    auto startTimer = timer::now();

    std::vector<std::vector<int>> waterToF(2048), expToF(2048);
    std::vector<int> tempVec(2048);
    for (int i = 0; i < 2048; i++)
        tempVec[i] = 0;
    for (int i = 0; i < 2048; i++) {
        waterToF[i] = tempVec;
        expToF[i] = tempVec;
    }

    std::ifstream waterFile(input_water);
    std::ifstream expFile(input_exp);
    for (int i = 0; i < 2048; i++) {
        for (int j = 0; j < 2048; j++) {
            waterFile >> waterToF[i][j];
            expFile >> expToF[i][j];
        }
    }

    detector detect(waterToF, expToF, info);
    detect.findingDifferences(6, 225);
    detect.printInFile(output);

    if (info) {
        std::chrono::duration<double> dur = timer::now() - startTimer;
        std::cout << "\t all time : " << dur.count() << std::endl;
    }
}

#endif //TIMEOFFLIGHT_ARGPARSE_H
