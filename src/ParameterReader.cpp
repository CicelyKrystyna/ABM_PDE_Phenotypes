#include "ParameterReader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

void ParameterReader::read(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Could not open file: " + filename);

    std::string line;
    std::string current_section;

    while (std::getline(infile, line)) {
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
        trim(line);
        if (line.empty()) continue;

        if (line.front() == '[' && line.back() == ']') {
            current_section = line.substr(1, line.size() - 2);
            trim(current_section);
            continue;
        }

        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        trim(key);
        trim(value);
        data[current_section][key] = value;
    }
}

std::string ParameterReader::get(const std::string& section, const std::string& key) const {
    if (!hasKey(section, key))
        throw std::runtime_error("Missing key [" + section + "] " + key);
    return data.at(section).at(key);
}

std::string ParameterReader::getString(const std::string& section, const std::string& key) const {
    return get(section, key);
}

std::string ParameterReader::getString(const std::string& section, const std::string& key, const std::string& defaultValue) const {
    return hasKey(section, key) ? data.at(section).at(key) : defaultValue;
}

int ParameterReader::getInt(const std::string& section, const std::string& key) const {
    return std::stoi(get(section, key));
}

int ParameterReader::getInt(const std::string& section, const std::string& key, int defaultValue) const {
    return hasKey(section, key) ? std::stoi(data.at(section).at(key)) : defaultValue;
}

double ParameterReader::getDouble(const std::string& section, const std::string& key) const {
    return std::stod(get(section, key));
}

double ParameterReader::getDouble(const std::string& section, const std::string& key, double defaultValue) const {
    return hasKey(section, key) ? std::stod(data.at(section).at(key)) : defaultValue;
}

std::vector<double> ParameterReader::getDoubleList(const std::string& section, const std::string& key, const std::vector<double>& defaults){
    if (!hasKey(section, key)) {
        return defaults;
    }

    std::string value = data.at(section).at(key); // raw string: "1.0, 2.5, 3.8"
    std::vector<double> result;
    std::stringstream ss(value);
    std::string item;

    while (std::getline(ss, item, ',')) {
        try {
            result.push_back(std::stod(item));
        } catch (...) {
            // fallback to 0.0 or throw exception if preferred
            result.push_back(0.0);
        }
    }

    return result.empty() ? defaults : result;
}


bool ParameterReader::hasKey(const std::string& section, const std::string& key) const {
    return data.count(section) && data.at(section).count(key);
}

void ParameterReader::trim(std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}
