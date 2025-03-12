#pragma once
#include <string>
#include <cstdint>


class Instruction {
public:
    Instruction() {};
    virtual ~Instruction() {};

    void run();
    virtual std::string to_string() = 0;

protected:

    virtual void decode() {};
    virtual void execute() {};

    void set_flag_Z(bool flag);
    void set_flag_Z(uint8_t value);
    void set_flag_Z(uint16_t value);

    void set_flag_N(bool flag);
    void set_flag_H(bool flag);
    void set_flag_C(bool flag);

    void set_flag_H_add(uint8_t x, uint8_t y, bool f=false); // bits 3 and 7
    void set_flag_H_add(uint16_t x, uint8_t y, bool f=false); // bit 3 and 7
    void set_flag_H_add(uint16_t x, int8_t y, bool f=false); // bit 3 and 7
    void set_flag_H_add(uint16_t x, uint16_t y, bool f=false); // bits 11 and 15
    void set_flag_C_add(uint8_t x, uint8_t y, bool f=false); // bits 3 and 7
    void set_flag_C_add(uint16_t x, uint8_t y, bool f=false); // bit 3 and 7
    void set_flag_C_add(uint16_t x, int8_t y, bool f=false); // bit 3 and 7
    void set_flag_C_add(uint16_t x, uint16_t y, bool f=false); // bits 11 and 15
    void set_flag_H_sub(uint8_t x, uint8_t y, bool f=false); // bits 3 and 7
    void set_flag_H_sub(uint16_t x, uint8_t y, bool f=false); // bit 3 and 7
    void set_flag_H_sub(uint16_t x, uint16_t y, bool f=false); // bits 11 and 15
    void set_flag_C_sub(uint8_t x, uint8_t y, bool f=false); // bits 3 and 7
    void set_flag_C_sub(uint16_t x, uint8_t y, bool f=false); // bit 3 and 7
    void set_flag_C_sub(uint16_t x, uint16_t y, bool f=false); // bits 11 and 15

};
