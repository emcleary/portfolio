#include "window_game.hpp"

#include "SDL2/SDL.h"

#include "../globals.hpp"


static const uint8_t NX = 20;
static const uint8_t NY = 18;
static const uint8_t NSPACE = 0;
static const uint8_t SCALE = 3;


WindowGame::WindowGame() : Window(NX, NY, NSPACE, SCALE) {}


void WindowGame::render() {
    SDL_memcpy(m_surface->pixels, m_buffer, 4 * m_buffer_size);
    m_pixel_index = 0;
    Window::render();
}


void WindowGame::set_pixel(uint8_t color) {
    if (m_pixel_index >= m_buffer_size) {
        std::cout << "Game screen pixel index out of bounds\n";
        std::cout << "index=" << +m_pixel_index << '\n';
        return;
    }

    m_buffer[m_pixel_index++] = IO::palette.m_palette[color];
}
