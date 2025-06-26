#include "ParameterReader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>

void ParameterReader::read(const std::string& filename) {
        std::ifstream infile(filename);
        if (!infile) throw std::runtime_error("Could not open file: " + filename);

        std::string line;
        std::string current_section;

        while (std::getline(infile, line)) {
            // Remove comments
            auto comment_pos = line.find('#');
            if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);

            // Trim
            trim(line);
            if (line.empty()) continue;

            // Section header
            if (line.front() == '[' && line.back() == ']') {
                current_section = line.substr(1, line.size() - 2);
                trim(current_section);
                continue;
            }

            // Key = value
            auto eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            std::string key = line.substr(0, eq_pos);
            std::string value = line.substr(eq_pos + 1);
            trim(key);
            trim(value);
            data[current_section][key] = value;
        }
    }

    std::string ParameterReader::getString(const std::string& section, const std::string& key) const {
        return get(section, key);
    }

    int ParameterReader::getInt(const std::string& section, const std::string& key) const {
        return std::stoi(get(section, key));
    }

    double ParameterReader::getDouble(const std::string& section, const std::string& key) const {
        return std::stod(get(section, key));
    }

    bool ParameterReader::hasKey(const std::string& section, const std::string& key) const {
        return data.count(section) && data.at(section).count(key);
    }

std::string ParameterReader::get(const std::string& section, const std::string& key) const {
        if (!hasKey(section, key))
            throw std::runtime_error("Missing key [" + section + "] " + key);
        return data.at(section).at(key);
    }

    

