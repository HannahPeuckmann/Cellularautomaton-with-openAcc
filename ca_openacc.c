// Cellular automaton parallelisation with OpenACC
// Hannah Peuckmann, Matr.-Nr.:791996
// Architectures and middleware for scientific computiong WiSe 21/22

/* (c) 1996,1997 Peter Sanders, Ingo Boesnach */
/* simulate a cellular automaton (serial version)
 * periodic boundaries
 *
 * #1: Number of lines
 * #2: Number of iterations to be simulated
 *
 */
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "random.h"
#include <fcntl.h>

/* horizontal size of the configuration */
#define XSIZE 1024

/* "ADT" State and line of states (plus border) */
typedef char State;
typedef State Line[XSIZE + 2];


#define TIME_GET(timer) \
	struct timespec timer; \
	clock_gettime(CLOCK_MONOTONIC, &timer)

#define TIME_DIFF(timer1, timer2) \
	((timer2.tv_sec * 1.0E+9 + timer2.tv_nsec) - \
	 (timer1.tv_sec * 1.0E+9 + timer1.tv_nsec)) / 1.0E+9

/* determine random integer between 0 and n-1 */
#define randInt(n) ((int)(nextRandomLEcuyer() * n))

/* --------------------- CA simulation -------------------------------- */

/* random starting configuration */
static void initConfig(Line *buf, int lines)
{
	int x, y;

	initRandomLEcuyer(424243);
	for (y = 1;  y <= lines;  y++) {
		for (x = 1;  x <= XSIZE;  x++) {
			buf[y][x] = randInt(100) >= 50;
		}
	}
}

/* a: pointer to array; x,y: coordinates; result: n-th element of anneal,
	    where n is the number of neighbors */
#define calc(a, x, y) \
    (a)[(y)-1][(x)-1] + (a)[(y)][(x)-1] + (a)[(y)+1][(x)-1] +\
		(a)[(y)-1][(x)  ] + (a)[(y)][(x)  ] + (a)[(y)+1][(x)  ] +\
		(a)[(y)-1][(x)+1] + (a)[(y)][(x)+1] + (a)[(y)+1][(x)+1]


/* treat torus like boundary conditions */
static void boundary_gpu(Line *buf, int lines)
{
	int x, y;

  #pragma acc parallel loop independent present(buf)
	for (y = 0;  y <= lines+1;  y++) {
		/* copy rightmost column to the buffer column 0 */
		buf[y][0      ] = buf[y][XSIZE];

		/* copy leftmost column to the buffer column XSIZE + 1 */
		buf[y][XSIZE+1] = buf[y][1    ];
	}

  #pragma acc parallel loop independent present(buf)
	for (x = 0;  x <= XSIZE+1;  x++) {
		/* copy bottommost row to buffer row 0 */
		buf[0][x      ] = buf[lines][x];

		/* copy topmost row to buffer row lines + 1 */
		buf[lines+1][x] = buf[1][x    ];
	}
}

/* make one simulation iteration with lines lines.
 * old configuration is in from, new one is written to to.
 */
static void simulate_gpu(Line *restrict from, Line *restrict to, int lines, State* anneal)
{
	int x, y;
  #pragma acc parallel loop present(from,to,anneal)
	for (y = 1;  y <= lines;  y++) {
		for (x = 1;  x <= XSIZE;  x++) {
      to[y][x]= anneal[calc(from, x, y)];
    }
	}
}

void write_matrix(Line* A, int lines, char file[]){
  int i, j;
  FILE *fp;

  fp = fopen(file, "w");

  for (i = 1; i <= lines; i++) {
      for (j = 1; j <= XSIZE; j++) {
          fprintf(fp, "%hhd ", A[i][j]);
      }
      fprintf(fp, "\n");
  }

  fclose(fp);
}

/* --------------------- measurement ---------------------------------- */

int main(int argc, char** argv)
{
	int lines, its, i;
	char* hash;
	Line *from, *to, *temp;

	/* init */

	assert(argc == 3);

	lines = atoi(argv[1]);
	its   = atoi(argv[2]);

	from = (Line*) calloc((lines + 2), sizeof(Line));
	to   = (Line*) calloc((lines + 2), sizeof(Line));

  if(to == NULL | from == NULL) {
    printf("ERROR: failure allocating device memory\n");
    exit(EXIT_FAILURE);
  }

	initConfig(from, lines);

	/* measurement loop */
	TIME_GET(start);

  State anneal[10] = {0, 0, 0, 0, 1, 0, 1, 1, 1, 1};
  const int matrixsize = (lines + 2) *  sizeof(Line);

  // copy 'from' to device and back out afterwards, create 'to' on the device
  #pragma acc data copy(from[0: matrixsize]), copyin(anneal[0:10]), create(to[0: matrixsize])
	for (i = 0;  i < its;  i++) {
    boundary_gpu(from, lines);
		simulate_gpu(from, to, lines, anneal);

	 	temp = from;  from = to;  to = temp;
	}
  #pragma end data

	TIME_GET(end);

	return EXIT_SUCCESS;
}