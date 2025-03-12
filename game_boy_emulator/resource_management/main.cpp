#include "SDL2/SDL.h"
#include "SDL2/SDL_ttf.h"

#include <memory>

#include "src/cartridge/cartridge.hpp"
#include "src/globals.hpp"


void usage() {
    std::cout << "usage:\n";
    std::cout << "   cassowary rom.gb\n";
}


int main(int argc, char **argv) {
    if (argc < 2) {
        usage();
        return -1;
    }

    std::string filename = argv[1];

    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);
    TTF_Init();

    std::shared_ptr<CartridgeInterface> cart = build_cart(filename);

    // initialize all global variables
    IO::state.initialize(filename, cart.get());
    CORE::cpu.initialize();
    SDL::ui.initialize();
    CORE::bus.initialize(cart);
    IO::timer.initialize(cart);
    IO::serial.initialize();
    IO::interrupts.initialize();
    CORE::registers.initialize();
    IO::lcds.initialize();
    IO::lcdc.initialize();
    IO::palette.initialize();
    IO::scroll.initialize();
    IO::dma.initialize();
    CORE::pipeline.initialize();
    CORE::ppu.initialize();
    SDL::gamepad.initialize();
    CORE::apu.initialize();

    while (!SDL::ui.get_die())
        CORE::cpu.step();

    SDL::ui.close();

    TTF_Quit();
    SDL_Quit();

    return 0;
}
