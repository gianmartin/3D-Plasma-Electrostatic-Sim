#ifndef _ARENA_H
#define _ARENA_H

#include <vector>
#include <cstddef>
#include <cassert>
#include <cstdint>

class ArenaAllocator {
public:
    ArenaAllocator(size_t size_bytes) : memory(size_bytes), offset(0) {}

    template <typename T>
    T* alloc(size_t count = 1) {
        size_t size_bytes = count * sizeof(T);
        size_t alignment = alignof(T);
        uintptr_t current_ptr = (uintptr_t)&memory[0] + offset;
        
        size_t padding = 0;
        size_t mask = alignment - 1;
        if (current_ptr & mask) {
            padding = alignment - (current_ptr & mask);
        }

        assert(offset + padding + size_bytes <= memory.size() && "Arena Out of Memory!");

        offset += padding;
        T* ptr = reinterpret_cast<T*>(&memory[offset]);
        offset += size_bytes;

        return ptr;
    }

    void reset() {offset = 0;}
    size_t get_used() const {return offset;}

private:
    std::vector<std::byte> memory;
    size_t offset;


};


#endif

