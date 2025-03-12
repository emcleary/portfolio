#pragma once

#include "arguments.hpp"
#include "instr_interface.hpp"


class CALL_N16 : public Instruction {
public:
    CALL_N16(BYTES::N16& src) : m_src(src) {}
        
    virtual std::string to_string() override;

protected:
    virtual void decode() final;
    virtual void execute() override;

    BYTES::N16& m_src;
};



class CALL_F_N16 final : public CALL_N16 {
public:
    CALL_F_N16(FLAG::FLAG& flag, BYTES::N16& src) : CALL_N16(src), m_flag(flag) {}

    virtual std::string to_string() final;

private:
    virtual void execute() final;
    
    FLAG::FLAG& m_flag;
};
