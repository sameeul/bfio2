#include <chrono>
#include <memory>
#include <vector>
#include <map>
#include <iostream>
//#include "ome_tiff_loader.h"
//#include "ome_zarr_loader.h"
#include "bioreader.h"
#include <sys/resource.h>
#include <pugixml.hpp>
#include <unistd.h>
#include <cstdlib>
#include <numeric>

// size_t partial_sum(std::shared_ptr<std::vector<tile_data_type>> data, size_t row_min, size_t row_max,  size_t col_min, size_t col_max, size_t col_width)
// {
//     size_t sum = 0;
//     for (int i=row_min; i<= row_max; ++i)
//     {
//         auto col_start_ptr = i*col_width;
//         for (int j=col_min; j<=col_max; j++)
//         {
//             sum = sum + data->at(col_start_ptr+j);
//         }
//     }
//     return sum;
// }



// void test1()
// {
//     std::cout<<"Test 1 - Virtual Tile From 4 Tiles" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     auto numRowTiles = imgLoader.GetRowTileCount();
//     auto numColTiles = imgLoader.GetColumnTileCount();  
//     auto start = std::chrono::steady_clock::now(); 

//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(0,0,0,0);
//     size_t sum = 0;
//     sum += partial_sum(tileData, 500, 1023,700,1023, 1024);
//     tileData = imgLoader.GetTileData(0,1,0,0);
//     sum += partial_sum(tileData, 500, 1023, 0,1070-1024, 56);
//     tileData = imgLoader.GetTileData(1,0,0);
//     sum += partial_sum(tileData, 0, 1050-1024,700,1023, 1024);
//     tileData = imgLoader.GetTileData(1,1,0,0);
//     sum += partial_sum(tileData, 0, 1050-1024,0,1070-1024,56);
    
//     std::cout << "Manual Total :" << sum << std::endl;

//     auto vTileData = imgLoader.GetVirtualTileData(Seq(500,1051), Seq(700, 1071), Seq(0,0), Seq(0,0) ,Seq(0,0));
//     size_t count = 0;
//     sum = 0;
//     for (auto x: *vTileData){
//         sum +=x;
//         count++;
//     }
//     std::cout <<"Virtual tile total: "<< sum<<std::endl;
//     auto end = std::chrono::steady_clock::now(); 
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

// }


// void test2()
// {
//     std::cout<<"Test 2 - Single Tile Subsection" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     auto numRowTiles = imgLoader.GetRowTileCount();
//     auto numColTiles = imgLoader.GetColumnTileCount();  
//     auto start = std::chrono::steady_clock::now(); 

//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(0,0,0);
//     size_t sum = 0;
//     sum += partial_sum(tileData, 0, 100,0,100, 1024);
   
//     std::cout << "Manual Total :" << sum << std::endl;

//     auto vTileData = imgLoader.GetVirtualTileData(Seq(0,101), Seq(0, 101), Seq(0,0), Seq(0,0) ,Seq(0,0));
//     size_t count = 0;
//     sum = 0;
//     for (auto x: *vTileData){
//         sum +=x;
//         count++;
//     }
//     std::cout <<"Virtual tile total: "<< sum<<std::endl;
//     auto end = std::chrono::steady_clock::now(); 
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

// }

// void test3()
// {
//     std::cout<<"Test 3 - Single Tile (Full)" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     auto numRowTiles = imgLoader.GetRowTileCount();
//     auto numColTiles = imgLoader.GetColumnTileCount();  
//     auto start = std::chrono::steady_clock::now(); 

//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(0,0,0);
//     size_t sum = 0;
//     sum += partial_sum(tileData, 0, 1023,0,1023, 1024);
   
//     std::cout << "Manual Total :" << sum << std::endl;

//     auto vTileData = imgLoader.GetVirtualTileData(Seq(0,1024), Seq(0, 1024), Seq(0,0), Seq(0,0) ,Seq(0,0));
//     size_t count = 0;
//     sum = 0;
//     for (auto x: *vTileData){
//         sum +=x;
//         count++;
//     }
//     std::cout <<"Virtual tile total: "<< sum<<std::endl;
//     auto end = std::chrono::steady_clock::now(); 
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

// }


// void test4()
// {
//     std::cout<<"Test 4 - Single Tile Memory usage" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     struct rusage rss1, rss2;
//     auto start = std::chrono::steady_clock::now(); 
//     auto tmp = getrusage(RUSAGE_SELF, &rss1);
//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(0,0,0);
//     tmp = getrusage(RUSAGE_SELF, &rss2);

//     std::cout<<"Memory usage for tile " << rss2.ru_maxrss - rss1.ru_maxrss << std::endl;
//     auto end = std::chrono::steady_clock::now(); 
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

// }

// // void test5()
// // {
// //     std::cout<<"Test 5 - Virtual Tile Stride" <<std::endl;
// //     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");

// //     std::shared_ptr<std::vector<tile_data_type>> tileData1 = imgLoader.GetBoundingBoxVirtualTileData(0,1079,0,1023,0,0);
// //     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetBoundingBoxVirtualTileDataStrideVersion(0,1079,1,0,1023,2,0,0,1);
// //     for(int i=0;i<10;i++){
// //         if (i%2==0){
// //             if(tileData1->at(i) != tileData->at(i/2)){std::cout<<"not match"<<std::endl;}
// //         }
// //     }
// //     size_t sum1 = std::accumulate(tileData->begin(), tileData->end(), 0);
// //     size_t sum2 = std::accumulate(tileData1->begin(), tileData1->end(), 0);
// //     std::cout << sum1 << " " << sum2<<std::endl;
// // }


// void test6()
// {
//     std::cout<<"Test 6 - Loop through tiles using index" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     for (int i=0; i<4; ++i){
//         std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileDataByIndex(i);
//         size_t sum = std::accumulate(tileData->begin(), tileData->end(), size_t(0));        
//         std::cout << "Tile id " << i << ", sum = "<< sum <<std::endl;
//     }
// }

// void test7()
// {
//     std::cout<<"Test 7 - Sanity check" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileDataByIndex(0);
// //    std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(0,0,0);

//     for (int i=0; i< imgLoader.GetTileHeight()*imgLoader.GetTileWidth(); ++i){
//         std::cout<<tileData->at(i)<<std::endl;
//     }

// }

void test8()
{
    std::cout<<"Test 8 - Iterator check" <<std::endl;
    BioReader<uint16_t>  imgLoader = BioReader<uint16_t> ("/home/ec2-user/data/3d-cell-viewer_ome_tiff_tiles.ome.tif",1);
    
    size_t tw = 1024;
    size_t th = 1024;
    size_t rs = 1024;
    size_t cs = 1024;
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

// void test9()
// {
//     std::cout<<"Test 6 - Loop through tiles using row and col" <<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/dev/imgloader/build/r01_x10_y05_z08.ome.tif");
//     for (int i=0; i<2; ++i){
//         for (int j=0; j<2; ++j){
//             std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetTileData(i,j,0);
//             size_t sum = std::accumulate(tileData->begin(), tileData->end(), size_t(0));        
//             std::cout << "Tile id " << i << ", sum = "<< sum <<std::endl;
//         }

//     }
// }

// void test10()
// {
//     std::cout<<"Test 10 - Validate CalcIFD"<<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/data/bfio_test_images/multi-channel-z-series.ome.tif");
//     std::cout<<imgLoader.CalcIFDIndex(2,0,0)<<std::endl;
//     std::shared_ptr<std::vector<tile_data_type>> tileData = imgLoader.GetVirtualTileData(Seq(0,166), Seq(0,438), Seq(1,2), Seq(0,0), Seq(0,0));
// }

// void test11(){
//     std::cout<<"Test 10 - Validate CalcIFD"<<std::endl;
//     BioReader imgLoader = BioReader("/mnt/hdd8/axle/data/polus-data/images/AICSImageIO/standard/tiff/actk_ome_tiff_tiles.ome.tif");
//     auto data = imgLoader.GetTileData(0,0,1,0,0);

// }

void test12(){
    std::cout<<"Test 12 - Read the whole image"<<std::endl;

    BioReader<float>  imgLoader = BioReader<float> ("/mnt/hdd8/axle/data/nyxus_demo_data/tiff/seg_dir/p01_x01_y01_wx0_wy0_c1.ome.tif",4);
    auto ih = imgLoader.GetImageHeight();
    auto iw = imgLoader.GetImageWidth();
    auto id = imgLoader.GetImageDepth();
    auto nc = imgLoader.GetChannelCount();
    auto nt = imgLoader.GetTstepCount();

    auto start = std::chrono::high_resolution_clock::now(); 
    std::shared_ptr<std::vector<float>> tileData = imgLoader.GetVirtualTileData(Seq(0,ih-1), Seq(0,iw-1), Seq(0,id-1), Seq(0,nc-1), Seq(0,nt-1));
    auto end = std::chrono::high_resolution_clock::now(); 
    size_t count = 0, sum = 0;
    for (auto x: *tileData){
        sum +=x;
        // count++;
        // std::cout<<x<<std::endl;
        // if (count==10) break;
    }
    std::cout <<"Virtual tile total: "<< sum<<std::endl;
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;
}

void test13(){

    std::shared_ptr<std::vector<uint16_t>> virtual_tile_data_1 = std::make_shared<std::vector<uint16_t>>(0) ;
	folly::resizeWithoutInitialization(*virtual_tile_data_1, 1024 * 1024 * 666);

    std::vector<uint16_t> src_data(1024*1024);
    for(auto i=0;i<1024*1024;i++){src_data[i] = rand()%100;}
    auto start = std::chrono::high_resolution_clock::now(); 
    for (auto i=0;i<666;++i){
        auto offset = i*1024*1024;
        memcpy((void *)(&(*(virtual_tile_data_1->data()+offset))), (void *)(&(src_data[0])), sizeof(uint16_t)*1024*1024);
    }
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

    std::shared_ptr<std::vector<uint16_t>> virtual_tile_data_2 = std::make_shared<std::vector<uint16_t>>(0) ;
	folly::resizeWithoutInitialization(*virtual_tile_data_2, 1024 * 1024 );
    start = std::chrono::high_resolution_clock::now(); 
        for (auto i=0;i<666;++i){
        auto offset = 0;
        memcpy((void *)(&(*(virtual_tile_data_2->data()+offset))), (void *)(&(src_data[0])), sizeof(uint16_t)*1024*1024);
    }
    end = std::chrono::high_resolution_clock::now(); 
    elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;

}


void test14(){
    std::cout<<"Test 14 - Read the whole image - OmeZarr"<<std::endl;

    BioReader<float>  imgLoader = BioReader<float> ("/mnt/hdd8/axle/data/nyxus_demo_data/zarr/seg_dir/p01_x01_y01_wx0_wy0_c1.ome.zarr",4);
    auto ih = imgLoader.GetImageHeight();
    auto iw = imgLoader.GetImageWidth();
    auto id = imgLoader.GetImageDepth();
    auto nc = imgLoader.GetChannelCount();
    auto nt = imgLoader.GetTstepCount();
    auto start = std::chrono::high_resolution_clock::now(); 
    std::shared_ptr<std::vector<float>> tileData = imgLoader.GetVirtualTileData(Seq(0,ih-1), Seq(0,iw-1), Seq(0,id-1), Seq(0,nc-1), Seq(0,nt-1));
    auto end = std::chrono::high_resolution_clock::now(); 
    size_t count = 0, sum = 0;
    for (auto x: *tileData){
        sum +=x;
    }
    std::cout <<"Virtual tile total: "<< sum<<std::endl;
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"elapsed_time " << elapsed_seconds.count() << std::endl;
}

int main(){
    // test1();
    // test2();
    // test3();
    // test4();
    // //test5();
    // test6();
    // //test7();
    // //test8();
    // test9();
    // test10();
    // test11();
    test12();
    //test13();
    test14();
    return 0;

}
