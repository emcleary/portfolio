#pragma once

#include "SDL2/SDL.h"

#include "../palette.hpp"


class Window {
public:
    Window(const uint8_t num_tiles_x, const uint8_t num_tiles_y,
            const uint8_t tile_space, const uint8_t scale)
            : m_num_tiles_x(num_tiles_x), m_num_tiles_y(num_tiles_y),
              m_tile_space(tile_space), m_scale(scale) {};

    ~Window() { close(); };

    void open();
    void close();
    void render();

    bool is_open() const { return m_is_running; }

    void set_window_title(std::string& title) { m_window_title = title; };
    void set_window_title(std::string&& title) { m_window_title = title; };

    SDL_Window* const get_window() const { return m_window; };
    SDL_Renderer* const get_renderer() const { return m_renderer; };
    SDL_Texture* const get_texture() const { return m_texture; };
    SDL_Surface* const get_surface() const { return m_surface; };

    const uint32_t get_window_width() const ;
    const uint32_t get_window_height() const ;
    const uint32_t get_surface_width() const ;
    const uint32_t get_surface_height() const ;

    const uint32_t get_tile_width() const { return 8; }
    const uint32_t get_tile_height() const { return 8; }
    const uint32_t get_tile_space() const { return m_tile_space; }

protected:

    uint32_t m_buffer_size;
    uint32_t* m_buffer;

    const uint8_t m_num_tiles_x;
    const uint8_t m_num_tiles_y;
    const uint8_t m_tile_space;
    uint32_t m_scale;

    bool m_is_running = false;

    std::string m_window_title;
    
    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    SDL_Texture* m_texture;
    SDL_Surface* m_surface;

    uint32_t m_tile_colors[4] = {
        COLOR::WHITE,
        COLOR::LIGHTGRAY,
        COLOR::DARKGRAY,
        COLOR::BLACK,
    };
};
