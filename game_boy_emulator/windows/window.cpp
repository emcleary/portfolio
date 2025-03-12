#include "window.hpp"

#include "SDL2/SDL.h"
#include "SDL2/SDL_ttf.h"

#include "../globals.hpp"
#include "../utilities.hpp"


void Window::open() {
    if (m_is_running)
        return;

    m_is_running = true;

    m_buffer_size = get_surface_height() * get_surface_width();
    m_buffer = new uint32_t[get_surface_height() * get_surface_width()];
    std::fill(m_buffer, m_buffer + m_buffer_size, 0xFF111111);

    SDL_CreateWindowAndRenderer(get_window_width(), get_window_height(), 0,
            &m_window, &m_renderer);

    m_surface = SDL_CreateRGBSurface(0, get_surface_width(), get_surface_height(), 32,
            0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);

    m_texture = SDL_CreateTexture(m_renderer, SDL_PIXELFORMAT_ARGB8888,
            SDL_TEXTUREACCESS_STREAMING,
            get_surface_width(), get_surface_height());

    SDL_SetWindowTitle(m_window, m_window_title.c_str());
}

void Window::close() {
    if (m_is_running) {
        SDL_FreeSurface(m_surface);
        SDL_DestroyTexture(m_texture);
        SDL_DestroyRenderer(m_renderer);
        SDL_DestroyWindow(m_window);
        delete[] m_buffer;

        m_surface = nullptr;
        m_texture = nullptr;
        m_renderer = nullptr;
        m_window = nullptr;
        m_is_running = false;
    }
}

void Window::render() {
    SDL_UpdateTexture(m_texture, NULL, m_surface->pixels, m_surface->pitch);
    SDL_RenderCopy(m_renderer, m_texture, NULL, NULL);
    SDL_RenderPresent(m_renderer);
}

const uint32_t Window::get_window_width() const {
    return m_scale * get_surface_width();
}

const uint32_t Window::get_window_height() const {
    return m_scale * get_surface_height();
}

const uint32_t Window::get_surface_width() const {
    uint32_t tile_width = get_tile_width();
    return m_num_tiles_x * tile_width + (m_num_tiles_x-1) * m_tile_space;
}

const uint32_t Window::get_surface_height() const {
    uint32_t tile_height = get_tile_height();
    return m_num_tiles_y * tile_height + (m_num_tiles_y-1) * m_tile_space;
}
