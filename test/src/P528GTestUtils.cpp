#include "P528GTest.h"
#include <fstream>
#include <sstream>
#include <filesystem>

void AppendDirectorySep(std::string &str) {
    #ifdef _WIN32
        str += "\\";
    #else
        str += "/";
    #endif
}

std::string GetDirectory(std::string name) {
    std::string dataDir(__FILE__);
    dataDir.resize(dataDir.find_last_of("/\\"));
    dataDir.resize(dataDir.find_last_of("/\\"));
    AppendDirectorySep(dataDir);
    dataDir += name;
    AppendDirectorySep(dataDir);
    return dataDir;
}

std::vector<InputsAndResult> ReadP528InputsAndResultFromPoints(const std::string& filename) {
    std::vector<InputsAndResult> testData;
    std::string dataDir = GetDirectory("data");
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

std::vector<std::string> GetDataTablesFileList()
{
    std::vector<std::string> ret;
    std::string dataTableDir = GetDirectory("Data Tables");
    for (const auto& entry : std::filesystem::directory_iterator(dataTableDir)) {
        //std::cout << entry.path() << std::endl;
        //std::string filename = entry.path().filename();
        ret.push_back(entry.path().filename().string());
    }
    return ret;
}

std::vector<InputsAndResult> ReadP528InputsAndResultFromDataTable(const std::string& filename, int testStep) {
    std::vector<InputsAndResult> testData;
    std::string dataDir = GetDirectory("Data Tables");
    std::ifstream file(dataDir + filename);
    InputsAndResult d;  // struct to store data from a single line of CSV

    double d__km;
    std::vector<double> h_1__meter;
    std::vector<double> h_2__meter;
    double f__mhz;
    int T_pol = 0;
    double p;

    //std::vector<double> A__db;
    double A_fs__db;

    std::vector<std::vector<std::string>> csvRows = readCSV(file);
    if (csvRows.size() <= 4) {
        return testData;
    }

    std::vector<std::string> cells;
   
    if (csvRows[0].size() > 0) {
        //std::cout << csvRows[0][0] << std::endl;
        string strCell = csvRows[0][0];
        strCell.resize(strCell.find("M"));
        f__mhz = std::stod(strCell);
        strCell = csvRows[0][0];
        strCell.resize(strCell.find(")"));
        strCell = strCell.substr(strCell.find("(") + 1);
        p = std::stod(strCell) * 100.0;
    }
    if (csvRows[1].size() > 2) {
        //std::cout << csvRows[1][2] << std::endl;
        for (int i = 2; i < csvRows[1].size(); i++) {
            h_2__meter.push_back(std::stod(csvRows[1][i]));
        }
    }
    if (csvRows[2].size() > 2) {
        //std::cout << csvRows[2][2] << std::endl;
        for (int i = 2; i < csvRows[2].size(); i++) {
            h_1__meter.push_back(std::stod(csvRows[2][i]));
        }
    }
    
    for (int r = 4; r < csvRows.size(); r += testStep) {
        //std::cout << csvRows[r][0] << std::endl;
        if (csvRows[r].size() > 2) {
            d__km = std::stod(csvRows[r][0]);
            A_fs__db = std::stod(csvRows[r][1]);

            for (int i = 2; i < csvRows[r].size(); i++) {
                d.d__km = d__km;
                d.h_1__meter = h_1__meter[i-2];
                d.h_2__meter = h_2__meter[i-2];
                d.f__mhz = f__mhz;
                d.T_pol = T_pol;
                d.p = p;

                if (i == 2) {
                    d.expectedResult.A_fs__db = A_fs__db;
                }
                else {
                    d.expectedResult.A_fs__db = NULL;
                }
                d.expectedResult.A__db = std::stod(csvRows[r][i]);
                testData.push_back(d);
            }
        }
    }

    return testData;
}

enum class CSVState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};

std::vector<std::string> readCSVRow(const std::string& row) {
    CSVState state = CSVState::UnquotedField;
    std::vector<std::string> fields{ "" };
    size_t i = 0; // index of the current field
    for (char c : row) {
        switch (state) {
        case CSVState::UnquotedField:
            switch (c) {
            case ',': // end of field
                fields.push_back(""); i++;
                break;
            case '"': state = CSVState::QuotedField;
                break;
            default:  fields[i].push_back(c);
                break;
            }
            break;
        case CSVState::QuotedField:
            switch (c) {
            case '"': state = CSVState::QuotedQuote;
                break;
            default:  fields[i].push_back(c);
                break;
            }
            break;
        case CSVState::QuotedQuote:
            switch (c) {
            case ',': // , after closing quote
                fields.push_back(""); i++;
                state = CSVState::UnquotedField;
                break;
            case '"': // "" -> "
                fields[i].push_back('"');
                state = CSVState::QuotedField;
                break;
            default:  // end of quote
                state = CSVState::UnquotedField;
                break;
            }
            break;
        }
    }
    return fields;
}

/// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
std::vector<std::vector<std::string>> readCSV(std::istream& in) {
    std::vector<std::vector<std::string>> table;
    std::string row;
    while (!in.eof()) {
        std::getline(in, row);
        if (in.bad() || in.fail()) {
            break;
        }
        auto fields = readCSVRow(row);
        table.push_back(fields);
    }
    return table;
}