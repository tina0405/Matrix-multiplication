#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>

#include <xmmintrin.h>

#define TEST_W 1024
#define TEST_H 1024

/* provide the implementations of naive_transpose,
 * sse_transpose, sse_prefetch_transpose
 */

#include "impl.c"
static long diff_in_us(struct timespec t1, struct timespec t2)
{
    struct timespec diff;
    if (t2.tv_nsec-t1.tv_nsec < 0) {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec - 1;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec + 1000000000;
    } else {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec;
    }
    return (diff.tv_sec * 1000000.0 + diff.tv_nsec / 1000.0);
}

int main()
{
   
    {
       

       int test_src1[TEST_W];
       int test_src2[TEST_H];
       for(int w=0;w<TEST_W;w++)
       {
           test_src1[w]=rand();
           test_src2[w]=rand();
       }
       int testout[TEST_H];
       struct timespec start, end;
       srand(time(NULL));



        clock_gettime(CLOCK_REALTIME, &start);
        naive_multiply(test_src1, test_src2, testout, 32, 32, 32, 32);
 	clock_gettime(CLOCK_REALTIME, &end);
        
        for (int y = 0; y <32 ; y++) {
            for (int x = 0; x < 32; x++)
                printf(" %2d", testout[y * 32 + x]);
            printf("\n");
        }
	printf("naive: \t\t %ld us\n", diff_in_us(start, end));

/*
#if defined(vertify) 
	 int test_src1[16] = { 0,  1,  2,  3,
                               4,  5,  6,  7,
                               8,  9, 10, 11,
                               12, 13, 14, 15
                            };
         int test_src2[16] = { 16, 17, 18, 19,
                               20, 21, 22, 23,
                               24, 25, 26, 27,
                               28, 29, 30, 31
                             };
         int testout[16];
         int expected[16] = { 152,  158,  164,  170,
                              504,  526,  548,  570,
                              856,  894,  932,  970,
                              1208, 1262, 1316, 1370
                            };       
        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++)
                printf(" %2d", test_src1[y * 4 + x]);
            printf("\n");
        }
        printf("\n");

        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++)
                printf(" %2d", test_src2[y * 4 + x]);
            printf("\n");
        }
        printf("\n");
        for (int y = 0; y <4 ; y++) {
            for (int x = 0; x < 4; x++)
                printf(" %2d", testout[y * 4 + x]);
            printf("\n");
        }

        assert(0 == memcmp(testout, expected, 16 * sizeof(int)) && "Verification fails");
#endif
*/
       
    }

    return 0;
}
