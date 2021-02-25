// Header only numeric CSV, based on Eigen
// Copyright (C) 2018  David Zuñiga-Noël <dzuniga at uma.es>

#ifndef CSV_H_
#define CSV_H_

// STL
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

// Eigen
#include <Eigen/Core>

namespace csv {

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> read(const std::string &path, char delim = ',') {

    std::string line;
    std::ifstream input(path);

    std::size_t cols = 0;
    std::vector<std::vector<Scalar>> data;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line.front() == '#') continue;

        std::vector<Scalar> data_row;
        data_row.reserve(cols);

        std::string value;
        std::stringstream line_stream(line);
        while (std::getline(line_stream, value, delim))
            data_row.push_back(static_cast<Scalar>(std::stod(value)));

        cols = std::max(cols, data_row.size());
        if (!data_row.empty()) data.push_back(data_row);
    }

    if (input.fail() && !input.eof())
        return Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(0, 0);

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix(data.size(), cols);
    for (std::size_t i = 0; i < data.size(); ++i) {
        std::vector<Scalar>& data_row = data.at(i);
        if (data_row.size() != cols) throw std::runtime_error("Data must be a table");
        matrix.row(i) = Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>>(data_row.data(), 1, cols); // TODO: is it safe? Does data get unallocated when exiting function (being data vectors cleared)?
    }

    return matrix;
}

template<typename Derived>
bool write(const Eigen::MatrixBase<Derived> &data, const std::string &path, int precision = 6, std::string delim = ",") {

    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, delim, "\n");

    std::ofstream output(path);

    output << std::setprecision(precision);
    output << std::fixed << data.format(fmt);

    return (!output.fail() && !output.bad());
}

} // namespace csv

#endif // CSV_H_

