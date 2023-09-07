#include "utilities.h"
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
#include <algorithm>
#include <regex>
#include "nlohmann/json.hpp"
// factory functions to create files, groups and datasets
#include "z5/factory.hxx"
// handles for z5 filesystem objects
#include "z5/filesystem/handle.hxx"
// attribute functionality
#include "z5/attributes.hxx"

bool CheckTileStatus(const std::string &fname){
	TIFF *tiff_ = TIFFOpen(fname.c_str(), "r");
	if (tiff_ != nullptr)
	{
		if (TIFFIsTiled(tiff_) == 0)
		{
			TIFFClose(tiff_);
			return false;
			} else
			{
			TIFFClose(tiff_);
			return true;
			}
	} else { throw (std::runtime_error("Tile Loader ERROR: The file can not be opened.")); }
}

std::string GetTiffType(const std::string &fname){
    TIFF * tiff = nullptr;
    tiff = TIFFOpen(fname.c_str(), "rm");
    std::string sample_type = "";
    if (tiff != nullptr){
        uint16_t sample_format = 0, bits_per_sample = 0, samples_per_pixel = 0;
        TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
        TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sample_format);
        // Interpret undefined data format as unsigned integer data
        if (sample_format < 1 || sample_format > 3) { sample_format = 1; }
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


std::string GetZarrType(const std::string &fname){
    std::unique_ptr<z5::filesystem::handle::File> zarr_ptr = std::make_unique<z5::filesystem::handle::File>(fname.c_str());
    nlohmann::json attributes;
    z5::readAttributes(*zarr_ptr, attributes);

    std::string metadata = attributes["metadata"].dump();
    std::regex type_regex("Type=\\\\\"(\\w+)");
    std::smatch matches;
    if(std::regex_search(metadata, matches, type_regex)) {
        if (matches[1].str() == "uint8") {return "uint8_t";}
        else if (matches[1].str() == "uint16") {return "uint16_t";}
        else if (matches[1].str() == "uint32") {return "uint32_t";}
        else if (matches[1].str() == "uint64") {return "uint64_t";}
        else if (matches[1].str() == "int8") {return "int8_t";}
        else if (matches[1].str() == "int16") {return "int16_t";}
        else if (matches[1].str() == "int32") {return "int32_t";}
        else if (matches[1].str() == "int64") {return "int64_t";}
        else if (matches[1].str() == "float") {return "float";}
        else if (matches[1].str() == "double") {return "double";}
        else {return "uint16_t";}
    }
    return "uint16_t";

}
