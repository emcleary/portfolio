
void Pipeline::process() {
    // A few lines of code omitted for clarity...
    (this->*m_method_fetcher)();
}

void Pipeline::fetcher_bgw_tile_start() {
    if (m_window_on) {
        if (IO::lcdc.get_win_enable()) {
            fetcher_w_tile_0();
            return;
        }
        m_window_on = false;
    }
    fetcher_bg_tile_0();
}

void Pipeline::fetcher_bg_tile_0() {
    m_tile_bg.set_map_addr();
    m_method_fetcher = &Pipeline::fetcher_bg_tile_1;
}

void Pipeline::fetcher_bg_tile_1() {
    m_method_fetcher = &Pipeline::fetcher_bg_data_lo_0;
}

void Pipeline::fetcher_bg_data_lo_0() {
    m_tile_bg.set_tile_index();
    m_tile_bg.set_tile_addr();
    m_tile_bg.load_data_lo();
    m_method_fetcher = &Pipeline::fetcher_bg_data_lo_1;
}

void Pipeline::fetcher_bg_data_lo_1() {
    m_method_fetcher = &Pipeline::fetcher_bg_data_hi_0;
}

void Pipeline::fetcher_bg_data_hi_0() {
    m_tile_bg.set_tile_index();
    m_tile_bg.set_tile_addr();
    m_tile_bg.load_data_hi();
    m_method_fetcher = &Pipeline::fetcher_bgw_data_hi_1;
}

void Pipeline::fetcher_bgw_data_hi_1() {
    if (can_fetch_obj() && (m_fifo_bgw.size() > 0)) {
        m_is_bgw_push_halted = true;
        fetcher_sprites_setup_0();
        return;
    }

    if (m_window_on_initial_fetch)
        m_window_on_initial_fetch = false;
    m_method_fetcher = &Pipeline::fetcher_bgw_push_to_fifo;
    fetcher_bgw_push_to_fifo();
}

void Pipeline::fetcher_bgw_push_to_fifo() {
    // only push to fifo if it is empty
    // otherwise come back on the next cycle and try again
    bool is_fifo_bgw_empty = m_fifo_bgw.size() == 0;
    if (is_fifo_bgw_empty) {
        for (int i = 0; i < 8; i++) {
            Pixel& pixel = m_window_on
                ? m_tile_w.get_pixel(i)
                : m_tile_bg.get_pixel(i);
            m_fifo_bgw.push(pixel);
        }

        if (m_window_on)
            m_tile_w.increment_tilemap_index_x();
        else
            m_tile_bg.increment_tilemap_index_x();

        m_method_fetcher = &Pipeline::fetcher_bgw_tile_start;
        if (m_window_on)
            return;
    }

    if (can_fetch_obj()) {
        m_is_bgw_push_halted = !is_fifo_bgw_empty;
        fetcher_sprites_setup_0();
        return;
    }
}
