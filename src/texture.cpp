#include "texture.h"
#include <framework/image.h>

#include <cmath>
#include <iostream>
#include <limits>



/*
    Called in "render.cpp"
    
    Method to compute the texels for the texture coordinates given for a set pixel. It receives a pointer to the image object 
    and its position in (i,j) format from texCoord. Since the center of the pixel starts at (0.5,0.5), after scaling the texel coordinates
    by the width and height of the image, we floor them to the nearest interger. The pixels are stored in an array and we compute their
    index similar to the PPM format back in class. 


*/
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    if (!features.enableTextureMapping) {
        return image.pixels[0];
    }

    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

   float width = image.width * texCoord.x;
   float height = image.height * texCoord.y;

   int index = floor(height) * image.width + floor(width);

    return image.pixels[index];
}


/*
    Bilinear Interpolation Computation of the Colors of a Texel

    The method works by taking the texel coordinate and computing the four texel nearest to it. That way
    it creates a square of the four pixels. From there the colors of the four pixels are all computed and
    the distances to each of the sides is turned into a weight. The closer to a pixel the texCoord is, the higher its weight.

    The index is computed in the same way as in the acquireTexel() method


*/
glm::vec3 acquireTexelBilinear(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    if (!features.enableTextureMapping || !features.extra.enableBilinearTextureFiltering) {
        return image.pixels[0];
    }


    float width = image.width * texCoord.x - 0.5f;
    float height = image.height * texCoord.y - 0.5f;


    //coordinates of square around point
    int A = floor(height) * image.width + floor(width);
    int B = floor(height) * image.width + ceil(width);
    int C = ceil(height) * image.width + floor(width);
    int D = ceil(height) * image.width + ceil(width);

    //colors of each of the four surrounding pixels 
    glm::vec3 cA = image.pixels[A];
    glm::vec3 cB = image.pixels[B];
    glm::vec3 cC = image.pixels[C];
    glm::vec3 cD = image.pixels[D];

    float wA = (ceil(width) - width) * (ceil(height) - height);
    float wB = (width - floor(width)) * (ceil(height) - height);
    float wC = (ceil(width) - width) * (height-floor(height));
    float wD = (width - floor(width)) * (height-floor(height));
    

    return wA * cA + wB * cB + wC * cC + wD * cD;
}



