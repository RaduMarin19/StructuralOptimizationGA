#pragma once

#include <vector>
#include <functional>

struct Chromosome {
	std::vector<bool> genes;

	bool operator==(const Chromosome& other) const {
		return genes == other.genes;
	}
};

namespace std {
	template <>
	struct hash<Chromosome> {
		size_t operator()(const Chromosome& c) const {
			size_t hash = 0;
			for (bool gene : c.genes) {
				hash <<= 1;
				hash ^= gene ? 1 : 0;
			}
			return hash;
		}
	};
}