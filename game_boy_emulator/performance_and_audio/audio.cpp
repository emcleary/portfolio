#include "audio.hpp"
#include <iostream>


static const int8_t MAX_SPEED = 10;
static const int8_t MIN_SPEED = -3;


Audio::~Audio() {
    SDL_CloseAudioDevice(device);
}

void Audio::initialize(uint32_t frequency, uint32_t sample_size, uint8_t channels) {
    m_default_frequency = frequency;

    SDL_zero(spec);
    spec.freq = frequency;
    spec.format = AUDIO_F32SYS;
    spec.channels = channels;
    spec.samples = sample_size;

    start_device();
}

void Audio::start_device() {
    SDL_CloseAudioDevice(device);
    device = SDL_OpenAudioDevice(nullptr, 0, &spec, nullptr, 0);
    SDL_PauseAudioDevice(device, false);
}

void Audio::queue(float* buffer, uint32_t size) {
    while (SDL_GetQueuedAudioSize(device)) {};
    SDL_QueueAudio(device, buffer, size * sizeof(float));
}

void Audio::slow_down() {
    if (m_speed_up_counter == MIN_SPEED) {
        std::cout << "Not allowed to go any slower!\n";
        return;
    }

    --m_speed_up_counter;
    scale_speed();
}

void Audio::speed_up() {
    if (m_speed_up_counter == MAX_SPEED) {
        std::cout << "Not allowed to go any faster!\n";
        return;
    }

    ++m_speed_up_counter;
    scale_speed();
}

void Audio::scale_speed() {
    int8_t scale = m_speed_up_counter >= 0
        ? m_speed_up_counter + 1
        : m_speed_up_counter - 1;

    spec.freq = scale > 0
        ? m_default_frequency * scale
        : m_default_frequency / (-scale);

    std::cout << "Speed changed to ";
    if (scale == 1)
        std::cout << "normal speed\n";
    else if (scale > 0)
        std::cout << +scale << "x speed\n";
    else
        std::cout << "1/" << +(-scale) << "x speed\n";

    start_device();
}
