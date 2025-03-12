#include "window_tiles.hpp"

#include "SDL2/SDL.h"

#include "../globals.hpp"


static const uint8_t NX = 16;
static const uint8_t NY = 24;
static const uint8_t NSPACE = 1;
static const uint8_t SCALE = 2;


WindowTiles::WindowTiles() : Window(NX, NY, NSPACE, SCALE) {};

WindowTiles::WindowTiles(const uint8_t num_tiles_x, const uint8_t num_tiles_y,
        const uint8_t tile_space, const uint8_t scale)
        : Window(num_tiles_x, num_tiles_y, tile_space, scale) {};


void WindowTiles::render() {
    if (!m_is_running)
        return;

    uint32_t tile_width = get_tile_width();
    uint32_t tile_height = get_tile_height();
    uint32_t tile_space = get_tile_space();
    uint16_t addr_tile = 0x8000;
    uint16_t x_draw = 0;
    uint16_t y_draw = 0;
    for (int y=0; y < m_num_tiles_y; y++) {
        for (int x=0; x < m_num_tiles_x; x++) {
            display_tile(addr_tile, x_draw, y_draw);
            addr_tile += 16;
            x_draw += tile_width + tile_space;
        }
        y_draw += tile_height + tile_space;
        x_draw = 0;
    }

    SDL_memcpy(m_surface->pixels, m_buffer, 4 * m_buffer_size);
    Window::render();
}


void WindowTiles::display_tile(uint16_t addr_start, uint16_t x, uint16_t y) {
    uint32_t line_width = get_surface_width();
    
    for (int tile_y = 0; tile_y < 16; tile_y += 2) {
        uint32_t index = x + (y + tile_y / 2) * line_width;
        uint16_t addr = addr_start + tile_y;
        uint8_t b0 = CORE::vram.read(addr);
        uint8_t b1 = CORE::vram.read(addr + 1);
        for (int bit = 7; bit >= 0; bit--) {
            uint8_t hi = (b0 >> bit) & 1;
            uint8_t lo = (b1 >> bit) & 1;
            uint8_t color = (hi << 1) | lo;
            m_buffer[index++] = m_tile_colors[color];
        }
    }
}
