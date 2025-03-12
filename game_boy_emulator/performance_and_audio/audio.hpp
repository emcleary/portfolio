#pragma once
#include "SDL2/SDL.h"


class Audio {
public:
    Audio() {};
    ~Audio();

    void initialize(uint32_t frequency, uint32_t sample_size, uint8_t channels);

    // Append queue AFTER queue empties --> controls frame rate
    void queue(float* buffer, uint32_t size);

    void slow_down();
    void speed_up();

private:
    void scale_speed();
    
    void start_device();

    SDL_AudioDeviceID device;
    SDL_AudioSpec spec;
    int8_t m_speed_up_counter = 0;
    uint32_t m_default_frequency;
};
