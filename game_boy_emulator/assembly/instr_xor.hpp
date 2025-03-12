#pragma once

#include "arguments.hpp"
#include "instr_interface.hpp"


template <typename T>
class XOR final : public Instruction {

public:
    XOR(REG::R8& dst, T& src) : m_dst(dst), m_src(src) {};
    virtual std::string to_string() final;

private:
    virtual void decode() final;
    virtual void execute() final;

    REG::R8& m_dst;
    T& m_src;
};


template class XOR<REG::R8>;
template class XOR<ADDR::R16>;
template class XOR<BYTES::N8>;

using XOR_R8_R8 = class XOR<REG::R8>;
using XOR_R8_AR16 = class XOR<ADDR::R16>;
using XOR_R8_N8 = class XOR<BYTES::N8>;
