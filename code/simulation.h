#pragma once

#include <fstream>
#include "colony.h"

void saveLifeHistory(const int generation, const Colony * const colonies[], const int sz, std::ofstream& file);
