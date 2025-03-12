#pragma once

#include <iostream>
#include <memory>

#include "cartridge/cartridge.hpp"
#include "memory.hpp"

class BUS {

public:
    BUS() {};
    ~BUS() {};

    void initialize(std::shared_ptr<Cartridge> cart);
    std::string to_string();

    uint8_t read(uint16_t addr);
    void write(uint16_t addr, uint8_t value);

    uint8_t read();
    uint16_t read16();
    uint16_t read16(uint16_t addr);
    void write16(uint16_t addr, uint16_t value);
    void increment(uint16_t addr); // included for incrementing DIV without resetting to 0

    void stack_push(uint8_t value);
    void stack_push16(uint16_t value);
    uint8_t stack_pop();
    uint16_t stack_pop16();

private:
    std::shared_ptr<Cartridge> m_cart;
};
