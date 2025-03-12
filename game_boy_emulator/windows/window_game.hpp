#pragma once

#include "window.hpp"


class WindowGame : public Window {
public:
    WindowGame();

    void render();
    void set_pixel(uint8_t color);

private:
    uint32_t m_pixel_index = 0;
};
