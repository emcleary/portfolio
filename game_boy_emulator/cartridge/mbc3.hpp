#pragma once

#include "cartridge_interface.hpp"


class MBC3 final : public Cartridge {
public:    
    MBC3() : Cartridge() {}

    virtual ~MBC3() {}

    virtual void initialize(Rom& rom) final;

    virtual uint8_t read(uint16_t addr) final;
    virtual void write(uint16_t addr, uint8_t value) final;

    virtual void tick() final;

    virtual void state_save(std::ofstream& file) final;
    virtual void state_load(std::ifstream& file) final;

private:

    virtual void load_battery() final;
    virtual void write_battery() final;

    void set_rom_bank_offsets();
    void set_ram_bank_index();

    bool m_ram_rtc_enable = false;
    uint8_t m_rom_bank_number = 1;
    uint8_t m_ram_bank_number = 0;
    uint8_t m_ram_bank_idx = 0;
    uint32_t m_rom_bank_offset = 0;
    bool m_latch_clock_data_written;

    RTC m_rtc, m_rtc_latched;
};
