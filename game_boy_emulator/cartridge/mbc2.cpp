#include "mbc2.hpp"


void MBC2::initialize(Rom& rom) {
    m_num_ram_banks = 1;
    Cartridge::initialize(rom);
    set_rom_bank_offsets();
}

uint8_t MBC2::read(uint16_t addr) {
    switch (addr & 0xE000) {
    case 0x0000:
    case 0x2000: // rom bank 0
        return m_rom[addr];
    case 0x4000:
    case 0x6000: // rom banks 0x01-0x0F
        return m_rom[(addr & 0x3FFF) | m_rom_bank_offset];
    case 0xA000:
        if (m_ram_enable) {
            uint8_t value = m_ram_banks[0].read(addr & 0xA1FF);
            return 0xF0 | value;
        }
        return 0xFF;
    default:
        std::cerr << std::format("Cartridge cannot read addr={:04X}\n", addr);
        exit(1);
    }
}

void MBC2::write(uint16_t addr, uint8_t value) {
    switch (addr & 0xE000) {
    case 0x0000:
    case 0x2000: {
        bool bit8set = (addr >> 8) & 0b1;
        if (bit8set) {
            m_rom_bank_number = value & 0x0F;
            if (m_rom_bank_number == 0)
                m_rom_bank_number = 1;
            set_rom_bank_offsets();
        } else {
            m_ram_enable = (value & 0b1111) == 0xA;
        }
        break;
    }
    case 0x4000:
    case 0x6000:
        break;
    case 0xA000:
        if (m_ram_enable) {
            m_ram_banks[0].write(addr & 0xA1FF, value & 0x0F);
            m_write_battery = true;
        }
        break;
    default:
        std::cerr << std::format("Cartridge cannot write to addr={:04X}\n", addr);
        exit(1);
    }
}

void MBC2::state_save(std::ofstream& file) {
    Cartridge::state_save(file);
    WRITE_BYTES(file, &m_ram_enable);
    WRITE_BYTES(file, &m_rom_bank_number);
    WRITE_BYTES(file, &m_rom_bank_offset);
}

void MBC2::state_load(std::ifstream& file) {
    Cartridge::state_load(file);
    READ_BYTES(file, &m_ram_enable);
    READ_BYTES(file, &m_rom_bank_number);
    READ_BYTES(file, &m_rom_bank_offset);
}

void MBC2::set_rom_bank_offsets() {
    m_rom_bank_offset  = m_rom_bank_number << 14;
    m_rom_bank_offset &= m_rom_size - 1;
}
