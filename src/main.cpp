// cpp standard library
#include <iostream>
#include <sstream>
#include <iomanip>
#include <array>
#include <limits>
#include <vector>
#include <algorithm>

// vendor library
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define HANDMADE_MATH_IMPLEMENTATION
#include "HandmadeMath.h"

#include "stable_fluids.hpp"

namespace
{
    const size_t image_width = 512;
    const size_t image_height = 512;
    // const float width_to_height_ratio = (float)image_width / (float)image_height;

    unsigned char image[3 * image_width * image_height];
}

static void write_to_image(size_t u, size_t v, hmm_vec3 value)
{
    value[0] = std::sqrt(value[0]);
    value[1] = std::sqrt(value[1]);
    value[2] = std::sqrt(value[2]);

    unsigned char r255 = static_cast<unsigned char>(255.99f * value[0]);
    unsigned char g255 = static_cast<unsigned char>(255.99f * value[1]);
    unsigned char b255 = static_cast<unsigned char>(255.99f * value[2]);

    unsigned char lower_bound = 0, upper_bound = 255;

    r255 = std::clamp(r255, lower_bound, upper_bound);
    g255 = std::clamp(g255, lower_bound, upper_bound);
    b255 = std::clamp(b255, lower_bound, upper_bound);

    image[v * 3 * image_width + u * 3 + 0] = r255;
    image[v * 3 * image_width + u * 3 + 1] = g255;
    image[v * 3 * image_width + u * 3 + 2] = b255;
}

static void output_image_to_file(std::string fileName)
{
    // std::cout << "fileName = " << fileName << std::endl;
    stbi_flip_vertically_on_write(true);
    if(stbi_write_png(fileName.c_str(), image_width, image_height, 3, image, 0) == 0)
    {
        std::cerr << "Error when saving image\n";
    }
}

int main(int argc, char* argv[])
{
    // for(size_t i = 0; i < image_width; ++i)
    // {
    //     for(size_t j = 0; j < image_height; ++j)
    //     {
    //         hmm_vec3 color{{ 0.0f, 0.0f, (float)(i + j) / (image_width + image_height) }};

    //         write_to_image(i, j, color);
    //     }
    // }
    // output_image_to_file("tmp.png");

    StableFluidSimulator simulator;
    simulator.Reset();

    for(size_t i = image_width/2 - 50; i < image_width/2 + 50; ++i)
    {
        for(size_t j = image_height/2 - 50; j < image_height/2 + 50; ++j)
        {
            // std::cout << i << ", " << j << std::endl;
            simulator.SetVX0(i, j, 15.0f * (2.0f * drand48() - 1.0f));
            simulator.SetVY0(i, j, 15.0f * (2.0f * drand48() - 1.0f));
            simulator.SetD0(i, j, 2.0f);
        }
    }
    simulator.AddSource();

    size_t iteration = 0;
    while(iteration++ < 1000)
    {
        simulator.CleanBuffer();
        simulator.VortConfinement();
        simulator.AnimateVel();
        simulator.AnimateDen();

        for(size_t i = 0; i < image_width; ++i)
        {
            for(size_t j = 0; j < image_height; ++j)
            {
                float d = simulator.GetBilinearFilteredDensity(i, j);

                hmm_vec3 color{{ 
                    1.0f - d, 
                    1.0f - d,
                    1.0f - d }};

                write_to_image(i, j, color);
            }
        }

        output_image_to_file("output/image" + std::to_string(iteration) + ".png");

    }

}