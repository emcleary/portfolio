#include "instr_xor.hpp"

#include <format>


template <typename T>
void XOR<T>::decode() {
    m_dst.decode();
    m_src.decode();
}

template <typename T>
void XOR<T>::execute() {
    uint8_t value_dst = m_dst.get();
    uint8_t value_src = m_src.get();
    uint8_t result = value_src ^ value_dst;
    m_dst.set(result);
    set_flag_Z(result == 0);
    set_flag_N(false);
    set_flag_H(false);
    set_flag_C(false);
}

template <typename T>
std::string XOR<T>::to_string() {
    return std::format("XOR {} {}", m_dst.to_string(), m_src.to_string());
}
