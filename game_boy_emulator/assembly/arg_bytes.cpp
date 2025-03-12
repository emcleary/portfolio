#include "arg_bytes.hpp"

#include "../globals.hpp"


namespace BYTES {


template <typename T>
T BYTES<T>::get() {
    return this->m_value;
}

template <typename T>
void BYTES<T>::set(T value) {
    std::cerr << "Bytes should never get set!\n";
    exit(1);
}


std::string N16::to_string() {
    return "N16";
}

void N16::decode() {
    m_value = CORE::bus.read16();
    IO::timer.add_cycles(2);
}


std::string N8::to_string() {
    return "N8";
}

void N8::decode() {
    m_value = CORE::bus.read();
    IO::timer.add_cycles(1);
}

std::string E8::to_string() {
    return "E8";
}

void E8::decode() {
    m_value = CORE::bus.read();
    IO::timer.add_cycles(1);
}


}
