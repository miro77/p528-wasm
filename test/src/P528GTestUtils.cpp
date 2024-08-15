#include "P528GTest.h"
#include <fstream>
#include <sstream>

void AppendDirectorySep(std::string &str) {
    #ifdef _WIN32
        str += "\\";
    #else
        str += "/";
    #endif
}

std::string GetDataDirectory() {
    std::string dataDir(__FILE__);
    dataDir.resize(dataDir.find_last_of("/\\"));
    dataDir.resize(dataDir.find_last_of("/\\"));
    AppendDirectorySep(dataDir);
    dataDir += "data";
    AppendDirectorySep(dataDir);
    return dataDir;
}

std::vector<InputsAndResult> ReadP528InputsAndResultTestData(const std::string& filename) {
    std::vector<InputsAndResult> testData;
    std::string dataDir = GetDataDirectory();
    std::ifstream file(dataDir + filename);
    std::string line;
    InputsAndResult d;  // struct to store data from a single line of CSV
    char c;  // single-character representing the comma (delimiter)
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (iss >> d.d__km >> c >> d.h_1__meter >> c >> d.h_2__meter >> c >> d.f__mhz >> c >> d.T_pol >> c >> d.p 
            >> c >> d.expectedReturn >> c >> d.expectedResult.propagation_mode >> c >> d.expectedResult.warnings >> c >> d.expectedResult.d__km
            >> c >> d.expectedResult.A__db >> c >> d.expectedResult.A_fs__db >> c >> d.expectedResult.A_a__db >> c >> d.expectedResult.theta_h1__rad
        ) {
            testData.push_back(d);
        }
    }
    return testData;
}
