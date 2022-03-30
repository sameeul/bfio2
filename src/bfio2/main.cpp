#include <chrono>
#include <memory>
#include <vector>
#include <map>
#include <iostream>
#include "ome_tiff_loader.h"
#include <sys/resource.h>
#include <pugixml.hpp>

#include <numeric>

size_t partial_sum(std::shared_ptr<std::vector<uint32_t>> data, size_t row_min, size_t row_max,  size_t col_min, size_t col_max, size_t col_width)
{
    size_t sum = 0;
    for (int i=row_min; i<= row_max; ++i)
    {
        auto col_start_ptr = i*col_width;
        for (int j=col_min; j<=col_max; j++)
        {
            sum = sum + data->at(col_start_ptr+j);
        }
    }
    return sum;
}



void test1()
{
    std::cout<<"Test 1 - Virtual Tile From 4 Tiles" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    auto numRowTiles = imgLoader.GetRowTileCount();
    auto numColTiles = imgLoader.GetColumnTileCount();  
    auto start = std::chrono::steady_clock::now(); 

    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(0,0);
    size_t sum = 0;
    sum += partial_sum(tileData, 500, 1023,700,1023, 1024);
    tileData = imgLoader.GetTileDataByRowCol(0,1);
    sum += partial_sum(tileData, 500, 1023, 0,1070-1024, 56);
    tileData = imgLoader.GetTileDataByRowCol(1,0);
    sum += partial_sum(tileData, 0, 1050-1024,700,1023, 1024);
    tileData = imgLoader.GetTileDataByRowCol(1,1);
    sum += partial_sum(tileData, 0, 1050-1024,0,1070-1024,56);
    
    std::cout << "Manual Total :" << sum << std::endl;

    auto vTileData = imgLoader.GetBoundingBoxVirtualTileData(500,1050, 700, 1070);
    size_t count = 0;
    sum = 0;
    for (auto x: *vTileData){
        sum +=x;
        count++;
    }
    std::cout <<"Virtual tile total: "<< sum<<std::endl;
    auto end = std::chrono::steady_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

}


void test2()
{
    std::cout<<"Test 2 - Single Tile Subsection" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    auto numRowTiles = imgLoader.GetRowTileCount();
    auto numColTiles = imgLoader.GetColumnTileCount();  
    auto start = std::chrono::steady_clock::now(); 

    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(0,0);
    size_t sum = 0;
    sum += partial_sum(tileData, 0, 100,0,100, 1024);
   
    std::cout << "Manual Total :" << sum << std::endl;

    auto vTileData = imgLoader.GetBoundingBoxVirtualTileData(0,100, 0, 100);
    size_t count = 0;
    sum = 0;
    for (auto x: *vTileData){
        sum +=x;
        count++;
    }
    std::cout <<"Virtual tile total: "<< sum<<std::endl;
    auto end = std::chrono::steady_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

}

void test3()
{
    std::cout<<"Test 3 - Single Tile (Full)" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    auto numRowTiles = imgLoader.GetRowTileCount();
    auto numColTiles = imgLoader.GetColumnTileCount();  
    auto start = std::chrono::steady_clock::now(); 

    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(0,0);
    size_t sum = 0;
    sum += partial_sum(tileData, 0, 1023,0,1023, 1024);
   
    std::cout << "Manual Total :" << sum << std::endl;

    auto vTileData = imgLoader.GetBoundingBoxVirtualTileData(0,1023, 0, 1023);
    size_t count = 0;
    sum = 0;
    for (auto x: *vTileData){
        sum +=x;
        count++;
    }
    std::cout <<"Virtual tile total: "<< sum<<std::endl;
    auto end = std::chrono::steady_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

}


void test4()
{
    std::cout<<"Test 4 - Single Tile Memory usage" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    struct rusage rss1, rss2;
    auto start = std::chrono::steady_clock::now(); 
    auto tmp = getrusage(RUSAGE_SELF, &rss1);
    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(0,0);
    tmp = getrusage(RUSAGE_SELF, &rss2);

    std::cout<<"Memory usage for tile " << rss2.ru_maxrss - rss1.ru_maxrss << std::endl;
    auto end = std::chrono::steady_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

}

void test5()
{
    std::cout<<"Test 5 - Virtual Tile Stride" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");

    std::shared_ptr<std::vector<uint32_t>> tileData1 = imgLoader.GetBoundingBoxVirtualTileData(0,1079,0,1023);
    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetBoundingBoxVirtualTileDataStrideVersion(0,1079,1,0,1023,2);
    for(int i=0;i<10;i++){
        if (i%2==0){
            if(tileData1->at(i) != tileData->at(i/2)){std::cout<<"not match"<<std::endl;}
        }
    }
    size_t sum1 = std::accumulate(tileData->begin(), tileData->end(), 0);
    size_t sum2 = std::accumulate(tileData1->begin(), tileData1->end(), 0);
    std::cout << sum1 << " " << sum2<<std::endl;
}


void test6()
{
    std::cout<<"Test 6 - Loop through tiles using index" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    for (int i=0; i<4; ++i){
        std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByIndex(i);
        size_t sum = std::accumulate(tileData->begin(), tileData->end(), size_t(0));        
        std::cout << "Tile id " << i << ", sum = "<< sum <<std::endl;
    }
}

void test7()
{
    std::cout<<"Test 7 - Sanity check" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByIndex(0);
//    std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(0,0);

    for (int i=0; i< imgLoader.GetTileHeight()*imgLoader.GetTileWidth(); ++i){
        std::cout<<tileData->at(i)<<std::endl;
    }

}

void test8()
{
    std::cout<<"Test 8 - Iterator check" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    
    size_t tw = 256;
    size_t th = 256;
    size_t rs = 200;
    size_t cs = 200;
    imgLoader.SetViewRequests(th, tw, rs, cs);
    auto ih = imgLoader.GetImageHeight();
    auto iw = imgLoader.GetImageWidth();
    for(size_t x=0; x<ih; x+=rs){
		size_t r_min = x;
		size_t r_max = x+rs-1;
		for(size_t y=0; y<iw; y+=cs){
			size_t c_min = y;
			size_t c_max = y+cs-1;
			auto tile_data = imgLoader.GetViewRequests(r_min, r_max, c_min, c_max);
            auto sum = std::accumulate(tile_data->begin(), tile_data->end(), size_t(0));
            std::cout << sum << std::endl;
		}
	}

}

void test9()
{
    std::cout<<"Test 6 - Loop through tiles using row and col" <<std::endl;
    OmeTiffLoader imgLoader = OmeTiffLoader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
    for (int i=0; i<2; ++i){
        for (int j=0; j<2; ++j){
            std::shared_ptr<std::vector<uint32_t>> tileData = imgLoader.GetTileDataByRowCol(i,j);
            size_t sum = std::accumulate(tileData->begin(), tileData->end(), size_t(0));        
            std::cout << "Tile id " << i << ", sum = "<< sum <<std::endl;
        }

    }
}
int main(){
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    //test7();
    //test8();
    test9();
    return 0;
}
