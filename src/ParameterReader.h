#ifndef PARAMETER_READER_H
#define PARAMETER_READER_H

#include <string>
#include <unordered_map>
#include <vector>

class ParameterReader {
public:
    void read(const std::string& filename);

    std::string getString(const std::string& section, const std::string& key) const;
    std::string getString(const std::string& section, const std::string& key, const std::string& defaultValue) const;

    int getInt(const std::string& section, const std::string& key) const;
    int getInt(const std::string& section, const std::string& key, int defaultValue) const;

    double getDouble(const std::string& section, const std::string& key) const;
    double getDouble(const std::string& section, const std::string& key, double defaultValue) const;
    std::vector<double> getDoubleList(const std::string& section, const std::string& key, const std::vector<double>& defaults);


    bool hasKey(const std::string& section, const std::string& key) const;

private:
    std::unordered_map<std::string, std::unordered_map<std::string, std::string> > data;

    std::string get(const std::string& section, const std::string& key) const;

    static void trim(std::string& s);
};

#endif // PARAMETER_READER_H
