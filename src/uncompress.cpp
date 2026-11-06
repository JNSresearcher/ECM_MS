#include <iostream>
#include <fstream> // Required for file stream operations
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm> // Required for std::transform
#include <zlib.h>
#include "./base64/base64.hpp"

    std::string decompress_zlib(const std::string& compressed_data) {
        z_stream zs;
        memset(&zs, 0, sizeof(zs));
        if (inflateInit(&zs) != Z_OK) {
            throw std::runtime_error("inflateInit failed while decompressing zlib data");
        }
        zs.next_in = (Bytef*)compressed_data.data();
        zs.avail_in = compressed_data.size();
        int ret;
        std::vector<char> uncompressed_data;
        do {
            char outbuffer[128];
            zs.next_out = (Bytef*)outbuffer;
            zs.avail_out = sizeof(outbuffer);
            ret = inflate(&zs, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR) {
                inflateEnd(&zs);
                throw std::runtime_error("inflate failed while decompressing zlib data");
            }
            uncompressed_data.insert(uncompressed_data.end(), outbuffer, outbuffer + sizeof(outbuffer) - zs.avail_out);
        } while (ret != Z_STREAM_END);
        inflateEnd(&zs);
        if (ret != Z_STREAM_END) {
            throw std::runtime_error("Unexpected end of stream while decompressing zlib data");
        }
        return std::string(uncompressed_data.begin(), uncompressed_data.end());
    }

extern "C" void c_routine (
                        char* input_text             // text to prepend to converted integer
                         )
{
    auto decoded_str = base64::from_base64(input_text);
    std::string decompressed_data = decompress_zlib(decoded_str);
    std::vector<int> intArray;
    intArray.clear();                                // Clear if using after Option 1
    intArray.resize(decompressed_data.length());     // Resize to avoid multiple allocations
    std::transform(decompressed_data.begin(), decompressed_data.end(), intArray.begin(),
                    [](char c) { return c - '0'; }); // Convert char '0' or '1' to int 0 or 1
    std::ofstream outputFile2("tmp.txt");            // Open file for writing
    if (outputFile2.is_open()) {
        for (int val : intArray) {
            outputFile2 << val+48 << " ";            // Write each integer followed by a space
        }
        outputFile2.close(); // Close the file
    } else {
        std::cerr << "Unable to open file tmp.txt for writing." << std::endl;
    }
}
