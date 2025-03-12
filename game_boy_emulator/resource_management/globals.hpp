#pragma once

#include "audio.hpp"
#include "apu.hpp"
#include "bus.hpp"
#include "cpu.hpp"
#include "dma.hpp"
#include "gamepad.hpp"
#include "interrupts.hpp"
#include "lcd_control.hpp"
#include "lcd_status.hpp"
#include "palette.hpp"
#include "pipeline.hpp"
#include "ppu.hpp"
#include "registers.hpp"
#include "scrolling.hpp"
#include "state.hpp"
#include "serial_data_transfer.hpp"
#include "timer.hpp"
#include "ui.hpp"


namespace SDL {
    inline Audio audio_player;
    inline Gamepad gamepad;
    inline UI ui;
}

namespace CORE {
    inline Memory vram = Memory(0x8000, 0x2000, 0);
    inline Memory wram = Memory(0xC000, 0x2000);
    inline Memory hram = Memory(0xFF80, 0x007F);
    inline Memory io = Memory(0xFEA0, 0x00E0);
    inline Memory oam = Memory(0xFE00, 0x00A0);

    inline BUS bus;
    inline CPU cpu;
    inline APU apu;
    inline PPU ppu;
    inline Pipeline pipeline;
    inline Registers registers;
}

namespace IO {
    inline DMA dma;
    inline LCD_Status lcds;
    inline LCD_Control lcdc;
    inline Timer timer;
    inline Serial serial;
    inline Interrupts interrupts;
    inline Scrolling scroll;
    inline Palette palette;
    inline State state;
}
