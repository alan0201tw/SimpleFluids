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
#undef STB_IMAGE_WRITE_IMPLEMENTATION

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#undef STB_IMAGE_IMPLEMENTATION

#define HANDMADE_MATH_IMPLEMENTATION
#include "HandmadeMath.h"
#undef HANDMADE_MATH_IMPLEMENTATION

#include "stable_fluids.hpp"

namespace
{
    const std::string texFileName = "resources/texture.png";
    const size_t image_width = 512;
    const size_t image_height = 512;
    
    unsigned char image[3 * image_width * image_height];

    int texWidth, texHeight;
    unsigned char* texData;
}

static hmm_vec3 sample_tex(size_t u, size_t v)
{
    return hmm_vec3{{
        (int)texData[3 * v * texWidth + 3 * u + 0] / 255.0f, 
        (int)texData[3 * v * texWidth + 3 * u + 1] / 255.0f,
        (int)texData[3 * v * texWidth + 3 * u + 2] / 255.0f }};
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
    int channel;
    texData = stbi_load(texFileName.c_str(), &texWidth, &texHeight, &channel, 0);
    std::cout << "SYS : Loading texture : " << texFileName 
        << ", width = " << texWidth << ", height = " << texHeight 
        << ", channel = " << channel << "\n";
    if(texData == nullptr)
    {
        std::cerr << "Error when loading texture : " << texFileName << "\n";
    }

    ////////////////////////////////

    StableFluidSimulator simulator(image_width, image_height);
    simulator.Reset();

    for(size_t i = image_width/2 - 50; i < image_width/2 + 50; ++i)
    {
        for(size_t j = image_height/2 - 50; j < image_height/2 + 50; ++j)
        {
            // std::cout << i << ", " << j << std::endl;
            simulator.SetVX0(i, j, 15.0f * (2.0f * drand48() - 1.0f));
            simulator.SetVY0(i, j, 15.0f * (2.0f * drand48() - 1.0f));
            // simulator.SetVX0(i, j, 1.0f);
            // simulator.SetVY0(i, j, 1.0f);
            simulator.SetD0(i, j, 10.0f);
        }
    }
    simulator.AddSource();

    size_t iteration = 0;
    
    while(iteration++ < 1000)
    {
        simulator.CleanBuffer();
        simulator.VortConfinement();
        simulator.AnimateVel();
        simulator.AnimateTex();
        simulator.AnimateDen();

        for(size_t i = 1; i < image_width-1; ++i)
        {
            for(size_t j = 1; j < image_height-1; ++j)
            {
                // float d = simulator.GetBilinearFilteredDensity(i, j);

                // hmm_vec3 color{{ 
                //     1.0f - d, 
                //     1.0f - d,
                //     1.0f - d }};

                auto uv = simulator.GetTextureCoord(i, j);
                // i, j is btw [0, image_width or image_height]

                size_t texCoordX = (uv.first / image_width) * (size_t)texWidth;
                size_t texCoordY = (uv.second / image_height) * (size_t)texHeight;

                // std::cout << texCoordX << ", " << texCoordY << "\n";

                hmm_vec3 color{{0.0f, 0.0f, 0.0f}};
                color += sample_tex(texCoordX, texCoordY);
                color += sample_tex(texCoordX+1, texCoordY);
                color += sample_tex(texCoordX, texCoordY+1);
                color += sample_tex(texCoordX+1, texCoordY+1);
                color *= 0.25f;

                // hmm_vec3 color{{
                //     idx / (512.0f * 512.0f * 3.0f), 
                //     idx / (512.0f * 512.0f * 3.0f),
                //     idx / (512.0f * 512.0f * 3.0f) }};

                write_to_image(i, j, color);
            }
        }

        output_image_to_file("output/image" + std::to_string(iteration) + ".png");
    }
}