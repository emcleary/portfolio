#include "factory.hpp"

#include "cartridge_types.hpp"
#include "header.hpp"
#include "mbc1.hpp"
#include "mbc2.hpp"
#include "mbc3.hpp"
#include "mbc5.hpp"
#include "nombc.hpp"
#include "rom.hpp"


std::shared_ptr<Cartridge> build_cart(std::string& filename) {
    Rom rom(filename);

    uint8_t num_ram_banks = rom.get_header().get_num_ram_banks();
    uint8_t num_rom_banks = rom.get_header().get_num_rom_banks();

    std::shared_ptr<Cartridge> cart = nullptr;
    
    switch(rom.get_header().cartridge_type) {
    case (CartTypes::NO_MBC):
        cart = std::make_shared<NoMBC>();
        break;
    case CartTypes::MBC1:
        cart = std::make_shared<MBC1>();
        break;
    case CartTypes::MBC1_RAM:
        cart = std::make_shared<MBC1>();
        cart->set_ram_banks(num_ram_banks);
        break;
    case CartTypes::MBC1_RAM_BATTERY:
        cart = std::make_shared<MBC1>();
        cart->set_ram_banks(num_ram_banks);
        cart->set_battery();
        break;
    case CartTypes::MBC2:
        cart = std::make_shared<MBC2>();
        break;
    case CartTypes::MBC2_BATTERY:
        cart = std::make_shared<MBC2>();
        cart->set_battery();
        break;
    case CartTypes::MBC3_TIMER_BATTERY:
        cart = std::make_shared<MBC3>();
        cart->set_timer();
        cart->set_battery();
        break;
    case CartTypes::MBC3_TIMER_RAM_BATTERY_2:
        cart = std::make_shared<MBC3>();
        cart->set_timer();
        cart->set_ram_banks(num_ram_banks);
        cart->set_battery();
        break;
    case CartTypes::MBC3:
        cart = std::make_shared<MBC3>();
        break;
    case CartTypes::MBC3_RAM_2:
        cart = std::make_shared<MBC3>();
        cart->set_ram_banks(num_ram_banks);
        break;
    case CartTypes::MBC3_RAM_BATTERY_2:
        cart = std::make_shared<MBC3>();
        cart->set_ram_banks(num_ram_banks);
        cart->set_battery();
        break;
    case CartTypes::MBC5:
        cart = std::make_shared<MBC5>();
        break;
    case CartTypes::MBC5_RAM:
        cart = std::make_shared<MBC5>();
        cart->set_ram_banks(num_ram_banks);
        break;
    case CartTypes::MBC5_RAM_BATTERY:
        cart = std::make_shared<MBC5>();
        cart->set_ram_banks(num_ram_banks);
        cart->set_battery();
        break;
    case CartTypes::MBC5_RUMBLE:
        cart = std::make_shared<MBC5>();
        break;
    case CartTypes::MBC5_RUMBLE_RAM:
        cart = std::make_shared<MBC5>();
        cart->set_ram_banks(num_ram_banks);
        break;
    case CartTypes::MBC5_RUMBLE_RAM_BATTERY:
        cart = std::make_shared<MBC5>();
        cart->set_ram_banks(num_ram_banks);
        cart->set_battery();
        cart->set_rumble();
        break;
    default:
        std::cerr << "Cartridge type " << rom.get_header().cartridge_type << " cannot be created\n";
        exit(1);
    }

    cart->set_rom_banks(num_rom_banks);
    cart->initialize(rom);
    return cart;
}
