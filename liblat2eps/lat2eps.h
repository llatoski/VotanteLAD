/*
 *  lat2eps 2.x
 *
 *  Copyright 2017 Andre R. de la Rocha
 *  
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
 *  
 *  Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES
 *  OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
 */

/** @file lat2eps.h */


#ifndef _LAT2EPS_H
#define _LAT2EPS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#define LAT2EPS_MAXQ    256   /*!< Maximum number of different site colors. */
#define LAT2EPS_MAXT   4096   /*!< Maximum number of text entries.          */
#define LAT2EPS_MAXL  16384   /*!< Maximum linear dimension of the lattice. */

#define LAT2EPS_VERS  "lat2eps 2.0"   /*!< Version string. */


/**
* Initializes the lattice resources. Must be called before any other lat2eps function.
* @param width  Lattice width (in sites).
* @param height Lattice height (in sites).
* @return       Zero for failure, non-zero for success.
*/
int lat2eps_init(unsigned int width, unsigned int height);


/**
* Releases the resources allocated by lat2eps_init().
*/
void lat2eps_release();


/**
* Sets the value of a lattice site.
* @param x Horizontal coordinate of the site.
* @param y Vertical coordinate of the site.
* @param s Value to set the site to.
*/
void lat2eps_set_site(unsigned int x, unsigned int y, int s);


/**
* Gets the value of a lattice site.
* @param x Horizontal coordinate of the site.
* @param y Vertical coordinate of the site.
* @return  Site value.
*/
int lat2eps_get_site(unsigned int x, unsigned int y);


/**
* Sets a color index to a palette entry defined in the 0xRRGGBB format.
* @param index Color index.
* @param pal   Palette entry in the 0xRRGGBB format.
*/
void lat2eps_set_color(unsigned int index, unsigned int pal);


/**
* Gets the palette definition associated with a color index.
* @param index Color index.
* @return      palette entry defined in the 0xRRGGBB format.
*/
unsigned int lat2eps_get_color(unsigned int index);


/**
* Adds a text message to the EPS output. Must be called before lat2eps_gen_eps().
* @param x      Horizontal coordinate where the text will be positioned. 0 is the leftmost coordinate, while the maximum value is defined by the lattice width.
* @param y      Vertical coordinate where the text will be positioned. 0 is the topmost coordinate, while the maximum value is defined by the lattice height.
* @param ax     Horizontal alignment. 0 for left-aligning the text relative to the X coordinate, 0.5 for centering it, 1 for right-aligning, etc.
* @param ay     Vertical alignment. 0 for placing the top of the text on the Y coordinate, 0.5 for centering it, 1 for placing the bottom of the text, etc.
* @param angle  Angle to rotate the text, in degrees (0 for horizontal left to right text).
* @param size   Font size.
* @param color  Color index.
* @param text   Text string.
* @note         Parentheses characters in the text must be escaped with backslashes.
*/
void lat2eps_text_out(float x, float y, float ax, float ay, float angle, unsigned int size, unsigned int color, const char *text);


/**
* Generates a lattice graphic in the EPS format.
* @param filename  Name of the EPS file that will be created, or NULL for outputting to stdout.
* @param xoff      Offset of the first lattice column that will be presented in the graphic.
* @param yoff      Offset of the first lattice row that will be presented in the graphic.
* @param width     Width (in sites) of the sublattice that will be presented in the graphic.
* @param height    Height (in sites) of the sublattice that will be presented in the graphic.
* @param border    Width of a black border that will be placed in the graphic (0 for no border).
* @param scale     Scale that will be used while generating the graphic (e.g., using 2 will create a 2x2 pixel square for each lattice site).
* @return          Zero for failure, non-zero for success.
*/
int lat2eps_gen_eps(const char *filename, unsigned int xoff, unsigned int yoff, unsigned int width, unsigned int height, unsigned int border, unsigned int scale);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAT2EPS_H */

