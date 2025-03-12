#pragma once

#include "window.hpp"


class WindowTiles : public Window {
public:
    WindowTiles();
    WindowTiles(const uint8_t num_tiles_x, const uint8_t num_tiles_y,
            const uint8_t tile_space, const uint8_t scale);
    
    void render();

    void display_tile(uint16_t add_start, uint16_t x, uint16_t y);

};
