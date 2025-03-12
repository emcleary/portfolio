#include "bus.hpp"

#include <format>

#include "addresses.hpp"
#include "globals.hpp"
// #include "utilities.hpp"


void BUS::initialize(std::shared_ptr<Cartridge> cart) {
    m_cart = cart;
}

std::string BUS::to_string() {
    uint16_t pc = CORE::registers.PC();
    return std::format("({:02X} {:02X} {:02X})", read(pc), read(pc+1), read(pc+2));
}

uint8_t BUS::read(uint16_t addr) {
    if (addr < 0x8000) 
        return m_cart->read(addr);

    else if ((0x8000 <= addr) && (addr < 0xA000)) {
        if (IO::lcds.is_mode_xfer())
            return 0xFF;
        return CORE::vram.read(addr);

    } else if ((0xA000 <= addr) && (addr < 0xC000))
        return m_cart->read(addr);

    else if ((0xC000 <= addr) && (addr < 0xE000))
        return CORE::wram.read(addr);

    else if ((0xE000 <= addr) && (addr < 0xFE00))
        return CORE::wram.read(addr & 0xDFFF);

    else if (0xFE00 <= addr && addr < 0xFEA0) {
        if (IO::dma.is_transferring()
            || IO::lcds.is_mode_oam()
            || IO::lcds.is_mode_xfer())
            return 0xFF;
        return CORE::oam.read(addr);

    } else if ((0xFF30 <= addr) && (addr < 0xFF40))
        return CORE::apu.channel3.read(addr);

    else if ((0xFF80 <= addr) && (addr < 0xFFFF))
        return CORE::hram.read(addr);

    switch (addr) {
    case ADDRESS::DIV:
        return IO::timer.get_div();
    case ADDRESS::TIMA:
        return IO::timer.get_tima();
    case ADDRESS::TMA:
        return IO::timer.get_tma();
    case ADDRESS::TAC:
        return IO::timer.get_tac();
    case ADDRESS::CONTROLLER:
        return SDL::gamepad.get_joyp();
    case ADDRESS::LCDC:
        return IO::lcdc.get_lcdc();
    case ADDRESS::LCDS:
        return IO::lcds.get_lcds();
    case ADDRESS::LY:
        return IO::lcds.get_ly();
    case ADDRESS::LYC:
        return IO::lcds.get_lyc();
    case ADDRESS::IF_REGISTER:
        return IO::interrupts.get_if();
    case ADDRESS::IE_REGISTER:
        return IO::interrupts.get_ie();
    case ADDRESS::BG_PALETTE:
        return IO::palette.get_bgp();
    case ADDRESS::OBJ_PALETTE_0:
        return IO::palette.get_obp0();
    case ADDRESS::OBJ_PALETTE_1:
        return IO::palette.get_obp1();
    case ADDRESS::SERIAL_DATA:
        return IO::serial.get_serial_transfer_data();
    case ADDRESS::SERIAL_CONTROL:
        return IO::serial.get_serial_transfer_control();
    case ADDRESS::SCROLL_X:
        return IO::scroll.get_scx();
    case ADDRESS::SCROLL_Y:
        return IO::scroll.get_scy();
    case ADDRESS::WINDOW_X_POS:
        return IO::scroll.get_win_x();
    case ADDRESS::WINDOW_Y_POS:
        return IO::scroll.get_win_y();
    case ADDRESS::NR10:
    case ADDRESS::NR11:
    case ADDRESS::NR12:
    case ADDRESS::NR13:
    case ADDRESS::NR14:
        return CORE::apu.channel1.read(addr);
    case ADDRESS::NR20:
    case ADDRESS::NR21:
    case ADDRESS::NR22:
    case ADDRESS::NR23:
    case ADDRESS::NR24:
        return CORE::apu.channel2.read(addr);
    case ADDRESS::NR30:
    case ADDRESS::NR31:
    case ADDRESS::NR32:
    case ADDRESS::NR33:
    case ADDRESS::NR34:
        return CORE::apu.channel3.read(addr);
    case ADDRESS::NR40:
    case ADDRESS::NR41:
    case ADDRESS::NR42:
    case ADDRESS::NR43:
    case ADDRESS::NR44:
        return CORE::apu.channel4.read(addr);
    case ADDRESS::NR50:
    case ADDRESS::NR51:
    case ADDRESS::NR52:
        return CORE::apu.read(addr);
    case 0xFF27:
    case 0xFF28:
    case 0xFF29:
    case 0xFF2A:
    case 0xFF2B:
    case 0xFF2C:
    case 0xFF2D:
    case 0xFF2E:
    case 0xFF2F:
        // unused memory
        return 0xFF;
    }

    // IO bytes not yet included above
    if ((0xFEA0 <= addr) && (addr < 0xFF80))
        return CORE::io.read(addr);

    std::cerr << std::format("Cannot read from addr={:04X}\n", addr);
    exit(1);
}

void BUS::write(uint16_t addr, uint8_t value) {
    if (addr < 0x8000) {
        m_cart->write(addr, value);
        return;

    } else if ((0x8000 <= addr) && (addr < 0xA000)) {
        if (!IO::lcds.is_mode_xfer())
            CORE::vram.write(addr, value);
        return;

    } else if ((0xA000 <= addr) && (addr < 0xC000)) {
        m_cart->write(addr, value);
        return;

    } else if ((0xC000 <= addr) && (addr < 0xE000)) {
        CORE::wram.write(addr, value);
        return;

    } else if ((0xE000 <= addr) && (addr < 0xFE00)) {
        CORE::wram.write(addr & 0xDFFF, value);
        return;

    } else if (0xFE00 <= addr && addr < 0xFEA0) {
        if (!IO::dma.is_transferring()
            && !IO::lcds.is_mode_oam()
            && !IO::lcds.is_mode_xfer())
            CORE::oam.write(addr, value);
        return;

    } else if ((0xFF30 <= addr) && (addr < 0xFF40)) {
        CORE::apu.channel3.write(addr, value);
        return;

    } else if ((0xFF80 <= addr) && (addr < 0xFFFF)) {
        CORE::hram.write(addr, value);
        return;
    }

    switch (addr) {
    case ADDRESS::DIV:
        IO::timer.set_div(value);
        return;
    case ADDRESS::TIMA:
        IO::timer.set_tima(value);
        return;
    case ADDRESS::TMA:
        IO::timer.set_tma(value);
        return;
    case ADDRESS::TAC:
        IO::timer.set_tac(value);
        return;
    case ADDRESS::DMA:
        IO::dma.start(value);
        break;
    case ADDRESS::CONTROLLER:
        SDL::gamepad.set_joyp(value);
        return;
    case ADDRESS::LCDC:
        IO::lcdc.set_lcdc(value);
        return;
    case ADDRESS::LCDS:
        IO::lcds.set_lcds(value);
        return;
    case ADDRESS::LY:
        IO::lcds.set_ly(value);
        return;
    case ADDRESS::LYC:
        IO::lcds.set_lyc(value);
        return;
    case ADDRESS::IF_REGISTER:
        IO::interrupts.set_if(value);
        return;
    case ADDRESS::IE_REGISTER:
        IO::interrupts.set_ie(value);
        return;
    case ADDRESS::BG_PALETTE:
        IO::palette.set_bgp(value);
        return;
    case ADDRESS::OBJ_PALETTE_0:
        IO::palette.set_obp0(value);
        return;
    case ADDRESS::OBJ_PALETTE_1:
        IO::palette.set_obp1(value);
        return;
    case ADDRESS::SERIAL_DATA:
        IO::serial.set_serial_transfer_data(value);
        return;
    case ADDRESS::SERIAL_CONTROL:
        IO::serial.set_serial_transfer_control(value);
        return;
    case ADDRESS::SCROLL_X:
        IO::scroll.set_scx(value);
        return;
    case ADDRESS::SCROLL_Y:
        IO::scroll.set_scy(value);
        return;
    case ADDRESS::WINDOW_X_POS:
        IO::scroll.set_win_x(value);
        return;
    case ADDRESS::WINDOW_Y_POS:
        IO::scroll.set_win_y(value);
        return;
    case ADDRESS::NR10:
    case ADDRESS::NR11:
    case ADDRESS::NR12:
    case ADDRESS::NR13:
    case ADDRESS::NR14:
        CORE::apu.channel1.write(addr, value);
        return;
    case ADDRESS::NR20:
    case ADDRESS::NR21:
    case ADDRESS::NR22:
    case ADDRESS::NR23:
    case ADDRESS::NR24:
        CORE::apu.channel2.write(addr, value);
        return;
    case ADDRESS::NR30:
    case ADDRESS::NR31:
    case ADDRESS::NR32:
    case ADDRESS::NR33:
    case ADDRESS::NR34:
        CORE::apu.channel3.write(addr, value);
        return;
    case ADDRESS::NR40:
    case ADDRESS::NR41:
    case ADDRESS::NR42:
    case ADDRESS::NR43:
    case ADDRESS::NR44:
        CORE::apu.channel4.write(addr, value);
        return;
    case ADDRESS::NR50:
    case ADDRESS::NR51:
    case ADDRESS::NR52:
        CORE::apu.write(addr, value);
        return;
    case 0xFF27:
    case 0xFF28:
    case 0xFF29:
    case 0xFF2A:
    case 0xFF2B:
    case 0xFF2C:
    case 0xFF2D:
    case 0xFF2E:
    case 0xFF2F:
        // unused memory
        return;
    }

    // IO bytes yet included above
    if ((0xFEA0 <= addr) && (addr < 0xFF80))
        return CORE::io.write(addr, value);

    std::cerr << std::format("Cannot write to addr={:04X}\n", addr);
    exit(1);
}

uint8_t BUS::read() {
    uint16_t addr = CORE::registers.increment_PC();
    return read(addr);
}

uint16_t BUS::read16() {
    uint8_t lo = read();
    uint8_t hi = read();
    return (hi << 8) | lo;
}

uint16_t BUS::read16(uint16_t addr) {
    uint8_t lo = read(addr);
    uint8_t hi = read(addr+1);
    return (hi << 8) | lo;
}

void BUS::increment(uint16_t addr) {
    uint8_t value = read(addr);
    value++;
    write(addr, value);
}

void BUS::write16(uint16_t addr, uint16_t value) {
    uint8_t lo = value & 0xFF;
    uint8_t hi = value >> 8;
    write(addr, lo);
    write(addr+1, hi);
}

void BUS::stack_push(uint8_t value) {
    CORE::registers.decrement_SP();
    uint16_t sp = CORE::registers.SP();
    write(sp, value);
}

void BUS::stack_push16(uint16_t value) {
    uint8_t hi = value >> 8;
    uint8_t lo = value & 0xFF;
    stack_push(hi);
    stack_push(lo);
}

uint8_t BUS::stack_pop() {
    uint16_t sp = CORE::registers.increment_SP();
    return read(sp);
}

uint16_t BUS::stack_pop16() {
    uint8_t lo = stack_pop();
    uint8_t hi = stack_pop();
    return (hi << 8) | lo;
}
