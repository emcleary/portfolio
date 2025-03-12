#pragma once

#include <memory>
#include <string>

#include "cartridge_interface.hpp"
#include "rom.hpp"

std::shared_ptr<Cartridge> build_cart(std::string& filename);
