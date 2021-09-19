#include <utility>
#include <cstdlib>

inline unsigned int hashInt(unsigned int a)
{
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) ^ (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a + 0xd3a2646c) ^ (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) ^ (a >> 16);
    return a;
}

inline unsigned hashIntPair(const std::pair<unsigned, unsigned> &p)
{
    unsigned h{hashInt(p.first)};
    h ^= hashInt(p.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
}

// For use in hash containers. For instance:
//   std::unordered_map<std::pair<int, int>, std::string, HashIntPairStruct>
//     int_to_string_map;
struct HashIntPairStruct
{
    size_t operator()(const std::pair<int, int> &p) const
    {
        return hashIntPair(p);
    }
};