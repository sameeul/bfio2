#include <string>
#ifdef __APPLE__
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#else
#include <tiffio.h>
#endif

std::string GetTiffType(const std::string &fname){
    TIFF * tiff = nullptr;
    tiff = TIFFOpen(fname.c_str(), "rm");
    std::string sample_type = "";
    if (tiff != nullptr){
        short sample_format = 0, bits_per_sample = 0, samples_per_pixel = 0;
        TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
        TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sample_format);
        TIFFClose(tiff);
        switch (sample_format) {
            case 1 :
                switch (bits_per_sample) {
                case 8: sample_type = "uint8_t";
                    break;
                case 16:sample_type = "uint16_t";
                    break;
                case 32:sample_type = "uint32_t";
                    break;
                case 64:sample_type = "uint64_t";
                    break;
                default:sample_type = "uint32_t";
                }
                break;
            case 2:
                switch (bits_per_sample) {
                case 8:sample_type = "int8_t";
                    break;
                case 16:sample_type = "int16_t";
                    break;
                case 32:sample_type = "int32_t";
                    break;
                case 64:sample_type = "int64_t";
                    break;
                default:sample_type = "uint32_t";
                }
                break;
            case 3:
                switch (bits_per_sample) {
                case 8:
                case 16:
                case 32:
                    sample_type = "float";
                    break;
                case 64:
                    sample_type = "double";
                    break;
                default:sample_type = "uint32_t";
                }
                break;
            default:sample_type = "uint32_t";
            }
    }

    return sample_type;
}