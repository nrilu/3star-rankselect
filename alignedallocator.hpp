//
// Created by nicco on 9/13/24.
//

#pragma once

#include <cstdint>
#include <new>
#include <limits>
#include <iostream>
#include <sys/mman.h>

template <typename T, std::size_t Alignment>
class AlignedAllocator {
public:
	using value_type = T;

	AlignedAllocator() noexcept = default;

	template <typename U>
		struct rebind {
		using other = AlignedAllocator<U, Alignment>;
	};

	bool operator==(const AlignedAllocator&) const noexcept {return true;} //needed for EXPL = EXPL_ to compile
	bool operator!=(const AlignedAllocator&) const noexcept {return false;}

	T* allocate(std::size_t n) {
		if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
			throw std::bad_alloc();

		size_t roundedSize = ((n * sizeof(T) + Alignment - 1) / Alignment) * Alignment;
		if (auto p = static_cast<T*>(aligned_alloc(Alignment, roundedSize))) {
			#if HUGEPAGES
				if(madvise(p, roundedSize, MADV_HUGEPAGE)!=0) {
					std::cerr << "madvise failed, but continuing anyway" << std::endl;
				}
			#endif
			return p;
		}
		throw std::bad_alloc();
	}

	void deallocate(T* p, std::size_t) noexcept {
		std::free(p);
	}
};

