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
#include "openssl/md5.h"
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

/* get MD5 checksum string of a memory chunk */
char* getMD5DigestStr(Line* buf, int lines)
{
	MD5_CTX ctx;
	unsigned char sum[MD5_DIGEST_LENGTH];
	int i;
	char* retval;
	char* ptr;

	MD5_Init(&ctx);
	for( i = 1; i<=lines; i++){	
		MD5_Update(&ctx, ((char*) buf[i])+1 , sizeof(State)*1024);
	}
	MD5_Final(sum, &ctx);

	retval = calloc(MD5_DIGEST_LENGTH * 2 + 1, sizeof(*retval));
	ptr = retval;

	for (i = 0; i < MD5_DIGEST_LENGTH; i++) {
		snprintf(ptr, 3, "%02X", sum[i]);
		ptr += 2;
	}

	return retval;
}

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

/* annealing rule from ChoDro96 page 34
 * the table is used to map the number of nonzero
 * states in the neighborhood to the new state
 */
static State anneal[10] = {0, 0, 0, 0, 1, 0, 1, 1, 1, 1};

/* a: pointer to array; x,y: coordinates; result: n-th element of anneal,
	    where n is the number of neighbors */
#define transition(a, x, y) \
	(anneal[(a)[(y)-1][(x)-1] + (a)[(y)][(x)-1] + (a)[(y)+1][(x)-1] +\
		(a)[(y)-1][(x)  ] + (a)[(y)][(x)  ] + (a)[(y)+1][(x)  ] +\
		(a)[(y)-1][(x)+1] + (a)[(y)][(x)+1] + (a)[(y)+1][(x)+1]])

/* treat torus like boundary conditions */
static void boundary_cpu(Line *buf, int lines)
{
	int x, y;
	for (y = 0;  y <= lines+1;  y++) {
		/* copy rightmost column to the buffer column 0 */
		buf[y][0      ] = buf[y][XSIZE];


		/* copy leftmost column to the buffer column XSIZE + 1 */
		buf[y][XSIZE+1] = buf[y][1    ];
	}

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
static void simulate_cpu(Line *from, Line *to, int lines)
{
	int x, y;

	boundary_cpu(from, lines);
	for (y = 1;  y <= lines;  y++) {
		for (x = 1;  x <= XSIZE;  x++) {
			to[y][x  ] = transition(from, x  , y);
		}
	}
}

/* treat torus like boundary conditions */
static void boundary_gpu(Line *buf, int lines)
{
	int x, y;
  #pragma acc data present(buf)
  #pragma acc parallel loop independent
	for (y = 0;  y <= lines+1;  y++) {
		/* copy rightmost column to the buffer column 0 */
		buf[y][0      ] = buf[y][XSIZE];

		/* copy leftmost column to the buffer column XSIZE + 1 */
		buf[y][XSIZE+1] = buf[y][1    ];
	}
  #pragma acc data present(buf)
  #pragma acc parallel loop independent
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
static void simulate_gpu(Line *from, Line *restrict to, int lines)
{
	int x, y;

  #pragma acc parallel loop present(from,to)
	for (y = 1;  y <= lines;  y++) {
    //#pragma acc loop independent
		for (x = 1;  x <= XSIZE;  x++) {
      static State anneal[10] = {0, 0, 0, 0, 1, 0, 1, 1, 1, 1};
      to[y][x]= anneal[(from)[(y)-1][(x)-1] + (from)[(y)][(x)-1] + (from)[(y)+1][(x)-1] +\
        (from)[(y)-1][(x)  ] + (from)[(y)][(x)  ] + (from)[(y)+1][(x)  ] +\
        (from)[(y)-1][(x)+1] + (from)[(y)][(x)+1] + (from)[(y)+1][(x)+1]];
		}
	}
}


/* --------------------- measurement ---------------------------------- */

int main(int argc, char** argv)
{
	int lines, its;
	int i;
	char* hash;
	Line *from, *to, *temp, *cputo, *cpufrom, *cputemp;

	/* init */

	assert(argc == 3);

	lines = atoi(argv[1]);
	its   = atoi(argv[2]);

	from = (Line*) calloc((lines + 2), sizeof(Line));
	to   = (Line*) calloc((lines + 2), sizeof(Line));
  cpufrom = (Line*) calloc((lines + 2), sizeof(Line));
	cputo   = (Line*) calloc((lines + 2), sizeof(Line));

	initConfig(from, lines);
  initConfig(cpufrom, lines);


	for (i = 0;  i < its;  i++) {
		simulate_cpu(cpufrom, cputo, lines);
      cputemp = cpufrom;
      cpufrom = cputo;
      cputo = cputemp;
	}
  hash = getMD5DigestStr(cpufrom, lines);
  printf("hash cpu: %s\n", hash);

	free(cpufrom);
	free(cputo);
	/* measurement loop */
	TIME_GET(start);
  // copy 'from' to device and back out afterwards, create 'to' on the device
  #pragma acc data copy(from), create(to)
	for (i = 0;  i < its;  i++) {
    boundary_gpu(from, lines);
		simulate_gpu(from, to, lines);
		temp = from;  from = to;  to = temp;
	} 

	TIME_GET(end);

	hash = getMD5DigestStr(from,lines);

	printf("hash gpu: %s\ttime: %.1f ms\n", hash, TIME_DIFF(start,end)*1000);

  free(from);
  free(to);
	free(hash);


	return EXIT_SUCCESS;
}
