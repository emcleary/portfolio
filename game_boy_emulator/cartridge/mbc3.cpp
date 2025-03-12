#include "mbc3.hpp"


void MBC3::initialize(Rom& rom) {
    Cartridge::initialize(rom);
    set_rom_bank_offsets();
    set_ram_bank_index();
}

uint8_t MBC3::read(uint16_t addr) {
    switch (addr & 0xE000) {
    case 0x0000:
    case 0x2000: // rom bank 0
        return m_rom[addr];
    case 0x4000:
    case 0x6000: // rom banks 0x01-0x7F
        return m_rom[(addr & 0x3FFF) | m_rom_bank_offset];
    case 0xA000: // ram
        if (m_ram_rtc_enable && (m_ram_bank_idx < 4))
            return m_ram_banks[m_ram_bank_idx].read(addr);
        else if (m_has_timer && (m_ram_bank_idx >= 8) && (m_ram_bank_idx <= 0xC)) {
            m_rtc_latched.set_registers(m_ram_bank_idx);
            return m_rtc_latched.read();
        }
        return 0xFF;
    default:
        std::cerr << std::format("Cartridge cannot read addr={:04X}\n", addr);
        exit(1);
    }
}

void MBC3::write(uint16_t addr, uint8_t value) {
    switch (addr & 0xE000) {
    case 0x0000:
        m_ram_rtc_enable = value == 0xA; // enables BOTH ram and RTC
        break;
    case 0x2000:
        m_rom_bank_number &= 0x7F;
        m_rom_bank_number = value == 0 ? 1 : value;
        set_rom_bank_offsets();
        break;
    case 0x4000:
        m_ram_bank_number = value;
        set_ram_bank_index();
        set_rom_bank_offsets();
        break;
    case 0x6000:
        if (m_latch_clock_data_written && (value == 1)) {
            if (m_has_timer)
                m_rtc_latched = m_rtc;
        }
        m_latch_clock_data_written = value == 0;
        break;
    case 0xA000:
        if (m_ram_rtc_enable) {
            if (m_ram_bank_idx < 4) {
                m_write_battery = true;
                m_ram_banks[m_ram_bank_idx].write(addr, value);
            } else if ((0x8 <= m_ram_bank_idx) && (m_ram_bank_idx <= 0xC)) {
                if (m_has_timer)
                    m_rtc_latched.write(value);
            }
        }
        break;
    default:
        std::cerr << std::format("Cartridge cannot write to addr={:04X}\n", addr);
        exit(1);
    }
}

void MBC3::tick() {
    Cartridge::tick();
    if (m_has_timer)
        m_rtc.increment();
}

void MBC3::state_save(std::ofstream& file)  {
    Cartridge::state_save(file);
    WRITE_BYTES(file, &m_ram_rtc_enable);
    WRITE_BYTES(file, &m_rom_bank_number);
    WRITE_BYTES(file, &m_ram_bank_number);
    WRITE_BYTES(file, &m_ram_bank_idx);
    WRITE_BYTES(file, &m_rom_bank_offset);
    WRITE_BYTES(file, &m_latch_clock_data_written);
    if (m_has_timer)
        file << m_rtc;
}

void MBC3::state_load(std::ifstream& file) {
    Cartridge::state_load(file);
    READ_BYTES(file, &m_ram_rtc_enable);
    READ_BYTES(file, &m_rom_bank_number);
    READ_BYTES(file, &m_ram_bank_number);
    READ_BYTES(file, &m_ram_bank_idx);
    READ_BYTES(file, &m_rom_bank_offset);
    READ_BYTES(file, &m_latch_clock_data_written);
    if (m_has_timer)
        file >> m_rtc;
}

void MBC3::load_battery() {
    if (m_has_battery) {
        std::cout << "Loading battery\n";
        m_battery->reset();
        for (auto& ram : m_ram_banks)
            *m_battery >> ram;
        if (m_has_timer)
            *m_battery >> m_rtc;
        m_battery->flush();
    }
}

void MBC3::write_battery() {
    if (m_has_battery && m_write_battery) {
        std::cout << "Writing battery\n";
        m_write_battery = false;
        m_battery->reset();
        for (auto& ram : m_ram_banks)
            *m_battery << ram;
        if (m_has_timer)
            *m_battery << m_rtc;
        m_battery->flush();
    }
}


void MBC3::set_rom_bank_offsets() {
    m_rom_bank_offset  = m_rom_bank_number << 14;
    m_rom_bank_offset &= m_rom_size - 1;
}

void MBC3::set_ram_bank_index() {
    m_ram_bank_idx = m_ram_bank_number;
}
