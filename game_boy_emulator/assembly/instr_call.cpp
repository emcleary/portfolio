#include "instr_call.hpp"

#include <format>

#include "../globals.hpp"


void CALL_N16::decode() {
    m_src.decode();
}

void CALL_N16::execute() {
    IO::timer.add_cycles(1);
    uint16_t pc = CORE::registers.PC();
    uint16_t value = m_src.get();
    CORE::bus.stack_push16(pc);
    IO::timer.add_cycles(2);
    CORE::registers.set_PC(value);
}

std::string CALL_N16::to_string() {
    return std::format("CALL {}", m_src.to_string());
}


void CALL_F_N16::execute() {
    if (m_flag.get())
        CALL_N16::execute();
}

std::string CALL_F_N16::to_string() {
    return std::format("CALL {},{}", m_flag.to_string(), m_src.to_string());
}
