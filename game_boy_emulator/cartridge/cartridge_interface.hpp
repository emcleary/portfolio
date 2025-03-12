#pragma once

#include <iostream>
#include <vector>

#include "battery.hpp"
#include "rom.hpp"

#include "../memory.hpp"


class Cartridge {
public:
    Cartridge() {}
    virtual ~Cartridge();

    void set_rom_banks(uint16_t n);
    void set_ram_banks(uint8_t n) { m_num_ram_banks = n; }
    void set_battery() { m_has_battery = true; }
    void set_timer() { m_has_timer = true; }
    void set_rumble() { m_has_rumble = true; }

    virtual void initialize(Rom& rom);

    virtual uint8_t read(uint16_t addr) = 0;
    virtual void write(uint16_t addr, uint8_t value) = 0;

    virtual void tick();

    virtual void state_save(std::ofstream& file);
    virtual void state_load(std::ifstream& file);

protected:

    virtual void load_battery();
    virtual void write_battery();

    uint16_t m_num_rom_banks = 1;
    uint8_t m_num_ram_banks = 0;
    uint8_t* m_rom = nullptr;
    size_t m_rom_size;
    bool m_has_battery = false;
    Battery* m_battery = nullptr;
    bool m_has_timer = false;
    bool m_has_rumble = false;
    std::vector<Memory> m_ram_banks;
    bool m_write_battery = false;
    uint32_t m_time = 0;
};
