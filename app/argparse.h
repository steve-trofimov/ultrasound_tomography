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

typedef std::chrono::high_resolution_clock timer;

//parametrs
const int size = 2048;
const int tickNumber = 3750;
const int numberOfWaves = 1;
const int signalDuration = 80;
const int delay = 80;
const int windowAkaike = 60;
const double cutoffLevel = 50;
bool timeTest = false;

namespace po = boost::program_options;

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
