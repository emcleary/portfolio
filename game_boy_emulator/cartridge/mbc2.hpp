#pragma once

#include "cartridge_interface.hpp"

#include <cassert>


// max 256 kB ROM, 512x4 bits RAM
class MBC2 final : public Cartridge {
public:    
    MBC2() : Cartridge() {};

    virtual ~MBC2() {};

    virtual void initialize(Rom& rom);

    virtual uint8_t read(uint16_t addr) final;
    virtual void write(uint16_t addr, uint8_t value) final;

    virtual void state_save(std::ofstream& file);
    virtual void state_load(std::ifstream& file);
    
private:
    void set_rom_bank_offsets();

    bool m_ram_enable = false;
    uint8_t m_rom_bank_number = 1;

    uint32_t m_rom_bank_offset = 0;
};
