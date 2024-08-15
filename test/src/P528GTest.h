#pragma once
#ifndef ITS_ITU_PSERIES_P528GTEST_H
#define ITS_ITU_PSERIES_P528GTEST_H

#include "gtest/gtest.h"

#include <p528.h>

struct InputsAndResult {
    double d__km;
    double h_1__meter;
    double h_2__meter;
    double f__mhz;
    int T_pol;
    double p;
    int expectedReturn;
    Result expectedResult;
};

std::vector<InputsAndResult> ReadP528InputsAndResultTestData(const std::string& filename);

#endif