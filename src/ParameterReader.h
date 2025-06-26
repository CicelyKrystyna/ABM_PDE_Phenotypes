#ifndef PARAMETER_READER_H
#define PARAMETER_READER_H

#include <string>
#include <unordered_map>

class ParameterReader {
public:
    // Load and parse a GetPot-style input file
    void read(const std::string& filename);

    // Accessors
    std::string getString(const std::string& section, const std::string& key) const;
    int getInt(const std::string& section, const std::string& key) const;
    double getDouble(const std::string& section, const std::string& key) const;

    // Check if a key exists in a section
    bool hasKey(const std::string& section, const std::string& key) const;

private:
    std::unordered_map<std::string, std::unordered_map<std::string, std::string> > data;

    // Internal value accessor with error checking
    std::string get(const std::string& section, const std::string& key) const;

    // Utility to trim whitespace from strings (in-place)
    static void trim(std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        size_t end = s.find_last_not_of(" \t\r\n");
        s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    };
};

#endif // PARAMETER_READER_H
