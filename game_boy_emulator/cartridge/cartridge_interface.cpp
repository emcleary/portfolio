#include "cartridge_interface.hpp"

#include <cstring>
#include <filesystem>


Cartridge::~Cartridge() {
    delete[] m_rom;
    if (m_battery) {
        write_battery();
        delete m_battery;
    }
}

void Cartridge::set_rom_banks(uint16_t n) {
    m_num_rom_banks = n;
    m_rom_size = 1024 * 16 * m_num_rom_banks;
}

void Cartridge::initialize(Rom& rom) {
    m_rom = new uint8_t[m_rom_size];
    memcpy(m_rom, rom.get_rom(), m_rom_size);

    for (int i=0; i < m_num_ram_banks; i++)
        m_ram_banks.push_back(Memory(0xA000, 0x2000));

    if (m_has_battery) {
        std::filesystem::path path(rom.get_filename());
        std::string battery_filename = path.replace_extension(".battery").string();
        
        m_battery = new Battery(battery_filename);
        if (!m_battery->new_save_file())
            load_battery();
    }
}

void Cartridge::tick() {
    if (++m_time % CPU_FREQUENCY == 0)
        write_battery();
}

void Cartridge::state_save(std::ofstream& file) {
    for (auto& bank : m_ram_banks)
        file << bank;
}

void Cartridge::state_load(std::ifstream& file) {
    for (auto& bank : m_ram_banks)
        file >> bank;
}

void Cartridge::load_battery() {
    if (m_has_battery) {
        std::cout << "Loading battery\n";
        m_battery->reset();
        for (auto& ram : m_ram_banks)
            *m_battery >> ram;
        m_battery->flush();
    }
}

void Cartridge::write_battery() {
    if (m_has_battery && m_write_battery) {
        std::cout << "Writing battery\n";
        m_write_battery = false;
        m_battery->reset();
        for (auto& ram : m_ram_banks)
            *m_battery << ram;
        m_battery->flush();
    }
}
