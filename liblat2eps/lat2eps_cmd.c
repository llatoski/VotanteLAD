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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "lat2eps.h"


#define LINE_BUFFER_LEN 16384
#define RC_FILE  ".lat2epsrc"


/* Shows usage info. */
void print_usage()
{
	fprintf(stderr, "Usage: lat2eps <xoff> <yoff> <width> <height> <border> <scale> [optional commands]  < input.dat  > output.eps\n");
}


/* Parses a buffer and returns pointers to tokens. */
unsigned int parse_buffer(char *buffer, unsigned int maxtokens, const char *separators, char **tokens, char **rest)
{
	unsigned int c = 0;

	if (buffer && separators && tokens) {
	
		char *tok = NULL;
		char *end = buffer + strlen(buffer);

		while ((c < maxtokens) && (tok = strsep(&buffer, separators))) {
			if (strlen(tok) > 0) {
				tokens[c++] = tok;
			}
		}

		if (rest) {
			if (tok && (tok + strlen(tok) + 1 < end)) {
				*rest = tok + strlen(tok) + 1;
			} else {
				*rest = NULL;
			}
		}
	}

	return c;
}


/* Processes an embedded command. */
void process_embedded_command(char *buffer)
{
	char *tokens[LAT2EPS_MAXQ + 10];
	const char *separators = ", \t\n\r";

	if (!strncasecmp(buffer, "TXT", 3) && strchr(separators, buffer[3])) {
		
		/* Text command */

		char *text = NULL, *end = NULL;
		unsigned int ntk = parse_buffer(buffer + 4, 7, separators, tokens, &text);

		if ((ntk == 7) && text) {

			float x = (float)atof(tokens[0]);
			float y = (float)atof(tokens[1]);
			float ax = (float)atof(tokens[2]);
			float ay = (float)atof(tokens[3]);
			float angle = (float)atof(tokens[4]);
			unsigned int size = (unsigned int)atoi(tokens[5]);
			unsigned int coloridx = (unsigned int)atoi(tokens[6]);
		
			while (isspace(*text)) ++text; /* trim left */
			end = text + strlen(text) - 1;
			while ((end >= text) && isspace(*end)) *end-- = 0;  /* trim right */

			lat2eps_text_out(x, y, ax, ay, angle, size, coloridx, text);
		}

	} else if (!strncasecmp(buffer, "COL", 3) && strchr(separators, buffer[3])) {
	
		/* Color command (changes a single color index) */

		unsigned int ntk = parse_buffer(buffer + 4, 2, separators, tokens, NULL);

		if (ntk == 2) {
			unsigned int coloridx = (unsigned int)atoi(tokens[0]);
			unsigned int pal = (unsigned int)strtol(tokens[1], NULL, 16);
			lat2eps_set_color(coloridx, pal);						
		}

	} else if (!strncasecmp(buffer, "PAL", 3) && strchr(separators, buffer[3])) {
	
		/* Palette command (changes the full palette) */

		unsigned int i;
		unsigned int ntk = parse_buffer(buffer + 4, LAT2EPS_MAXQ, separators, tokens, NULL);
	
		for (i = 0; i < ntk; ++i) {
			unsigned int pal = (unsigned int)strtol(tokens[i], NULL, 16);
			lat2eps_set_color(i, pal);						
		}
	}
}


/* Reads the lattice data. */
void read_lattice(FILE *f, int xoff, int yoff, int width, int height)
{
	char buffer[LINE_BUFFER_LEN];
	
	while (fgets(buffer, sizeof(buffer), f)) {

		char *p = buffer;
		while (isspace(*p)) ++p;

		if (*p) {

			if (*p == '#') {    /* Lines beginning with # are comments that may contain embedded commands. */

				process_embedded_command(p + 1);

			} else {    /* Handle site data */

				char *tokens[3];
				unsigned int ntk = parse_buffer(p, 3, " \t\n\r", tokens, NULL);

				if (ntk == 3) {
					int x = atoi(tokens[0]);
					int y = atoi(tokens[1]);
					int s = atoi(tokens[2]);
					if ((x >= xoff) && (x < xoff + width) && (y >= yoff) && (y < yoff + height)) {
						lat2eps_set_site(x-xoff, y-yoff, s);
					}
				}
			}
		}
	}
}


/* Reads default configuration defined through embedded commands in the ~/.lat2epsrc file. */
void read_rc_defaults()
{
	char buffer[LINE_BUFFER_LEN];
	char filename[PATH_MAX];
	char *home;
	
	if ((home = getenv("HOME"))) {

		FILE *f;
		snprintf(filename, sizeof(filename), "%s/%s", home, RC_FILE);
		
		if ((f = fopen(filename, "r"))) {

			while (fgets(buffer, sizeof(buffer), f)) {

				char *p = buffer;
				while (isspace(*p)) ++p;

				if (*p == '#') {
					/* Lines beginning with # are comments that may contain embedded commands. */
					process_embedded_command(p + 1);
				}
			}
		
			fclose(f);
		}
	}
}


int main(int argc, char *argv[])
{
	int argidx, xoff, yoff, width, height, border, scale;

	if (argc < 7) {
		print_usage();
		return -1;
	}
	
	argidx = 1;
	
	/* Gets required arguments */
	xoff = atoi(argv[argidx++]);
	yoff = atoi(argv[argidx++]);
	width = atoi(argv[argidx++]);
	height = atoi(argv[argidx++]);
	border = atoi(argv[argidx++]);
	scale = atoi(argv[argidx++]);
	
	if ((width <= 0) || (width > LAT2EPS_MAXL) || (height <= 0) || (height > LAT2EPS_MAXL) || (border < 0) || (scale <= 0)) {
		print_usage();
		return -1;
	}

	/* Initializes lat2eps lib */
	if (!lat2eps_init(width, height)) {
		fprintf(stderr, "Initialization error.\n");
		return -1;
	}

	read_rc_defaults();  /* Reads configuration (embedded commands) from rc file. */

	read_lattice(stdin, xoff, yoff, width, height);  /* Reads lattice data (and embedded commands) from stdin. */

	/* Treats remaining arguments after required ones as possible embedded commands. */
	while (argidx < argc) {
		process_embedded_command(argv[argidx++]);
	}

	/* Outputs EPS data to stdout */
	if (!lat2eps_gen_eps(NULL, 0, 0, width, height, border, scale)) {
		lat2eps_release();
		fprintf(stderr, "EPS generation error.\n");
		return -1;
	}

	/* Releases resources allocated by lat2eps lib. */
	lat2eps_release();

	return 0;
}

