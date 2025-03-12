#pragma once

#include <cinttypes>
#include <fstream>
#include <string>

#include "rtc.hpp"
#include "../memory.hpp"


class Battery {
public:
    Battery(std::string& filename);
    Battery(std::string&& filename);

    ~Battery();

    void reset();
    void flush();

    // TODO: make as template?
    void operator<<(Memory& memory);
    void operator>>(Memory& memory);
    void operator<<(RTC& rtc);
    void operator>>(RTC& rtc);
    
    bool new_save_file() const;

private:
    void initialize();
    
    std::string m_filename;
    std::fstream m_file;
    bool m_new_file = false;
};
