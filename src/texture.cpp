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