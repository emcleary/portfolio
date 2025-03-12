#include "window_map.hpp"

#include "SDL2/SDL.h"

#include "../globals.hpp"


static const uint8_t NX = 32;
static const uint8_t NY = 64;
static const uint8_t NSPACE = 1;
static const uint8_t SCALE = 1;


WindowMap::WindowMap() : WindowTiles(NX, NY, NSPACE, SCALE) {}


void WindowMap::render() {
    if (!m_is_running)
        return;

    uint32_t tile_width = get_tile_width();
    uint32_t tile_height = get_tile_height();
    uint32_t tile_space = get_tile_space();
    uint16_t addr = 0x9800;
    uint16_t addr_tile = IO::lcdc.get_bgw_data_area();
    uint8_t tile_offset = addr_tile == 0x8000 ? 0 : 0x80;
    uint16_t x_draw = 0;
    uint16_t y_draw = 0;
    for (int y=0; y < m_num_tiles_y; y++) {
        for (int x=0; x < m_num_tiles_x; x++) {
            uint16_t tile_num = CORE::vram.read(addr++) + tile_offset;
            display_tile(addr_tile + 16*tile_num, x_draw, y_draw);
            x_draw += tile_width + tile_space;
        }
        y_draw += tile_height + tile_space;
        x_draw = 0;
    }

    SDL_memcpy(m_surface->pixels, m_buffer, 4 * m_buffer_size);
    Window::render();
}
