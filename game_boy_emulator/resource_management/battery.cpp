#include "battery.hpp"


Battery::Battery(std::string& filename) : m_filename(filename) {
    initialize();
}

Battery::Battery(std::string&& filename) : m_filename(filename) {
    initialize();
}

Battery::~Battery() {
    if (m_file.is_open())
        m_file.close();
}

void Battery::operator<<(Memory& memory) {
    m_file << memory;
}

void Battery::operator>>(Memory& memory) {
    m_file >> memory;
}

void Battery::operator<<(RTC& rtc) {
    m_file << rtc;
}

void Battery::operator>>(RTC& rtc) {
    m_file >> rtc;
}

void Battery::reset() {
    m_file.seekp(0, std::fstream::beg);
}

void Battery::flush() {
    m_file.flush();
}
    
bool Battery::new_save_file() const {
    return m_new_file;
}

void Battery::initialize() {
    std::fstream init(m_filename, std::fstream::binary | std::fstream::app);
    if (init.is_open()) {
        if (init.tellp() == 0) {
            std::cout << "Created new battery file " << m_filename << '\n';
            m_new_file = true;
        } else {
            std::cout << "Using existing battery file " << m_filename << '\n';
            m_new_file = false;
        }
    } else {
        std::cerr << "ERROR creating/loading battery file " << m_filename << '\n';
        exit(-1);
    }
    init.close();

    m_file.open(m_filename, std::fstream::binary | std::fstream::in | std::fstream::out);
}
