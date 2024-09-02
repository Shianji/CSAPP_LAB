/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]) {
  int i, j, ii, jj, a1, a2, a3, a4, a5, a6, a7, a0;
  if (M == 32) {

    for (i = 0; i < N; i += 8) {
      for (j = 0; j < M; j += 8) {
        for (ii = i; ii < i + 8; ii++) {
          jj = j;
          a0 = A[ii][jj];
          a1 = A[ii][jj + 1];
          a2 = A[ii][jj + 2];
          a3 = A[ii][jj + 3];
          a4 = A[ii][jj + 4];
          a5 = A[ii][jj + 5];
          a6 = A[ii][jj + 6];
          a7 = A[ii][jj + 7];
          B[jj][ii] = a0;
          B[jj + 1][ii] = a1;
          B[jj + 2][ii] = a2;
          B[jj + 3][ii] = a3;
          B[jj + 4][ii] = a4;
          B[jj + 5][ii] = a5;
          B[jj + 6][ii] = a6;
          B[jj + 7][ii] = a7;
        }
      }
    }
  } else if (M == 64) {
    int i, j;

    for (i = 0; i < N; i += 8) {      // 枚举每八行
        for (j = 0; j < M; j += 8) {  // 枚举每八列
            // 这里用这些临时变量，如果你查看过A和B的地址的话，你会发现A和B的地址差距是64的整数倍（0x40000），
            // 那么 直接赋值的话，在对角线的时候 每一个Load A[i][i]紧跟Store B[i][i],将造成比较多的
            // eviction
            int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, cnt;
            // 1块8*8，我们分成4块来做，每一块4*4
            // 这是 左上，且将左下的块移动到 B中的右上，这样是为了更高效地利用cache
            for (cnt = 0; cnt < 4; ++cnt, ++i) {  // 枚举0~8中的每一行，一行八列
                temp1 = A[i][j];  // 这样我们就一次取出了8个元素，我们的A的miss就只有一次了 原始4*4
                                  // 则是两次
                temp2 = A[i][j + 1];
                temp3 = A[i][j + 2];
                temp4 = A[i][j + 3];  // 左上

                temp5 = A[i][j + 4];  // 右上
                temp6 = A[i][j + 5];
                temp7 = A[i][j + 6];
                temp8 = A[i][j + 7];

                B[j][i] = temp1;  // 左上翻转
                B[j + 1][i] = temp2;
                B[j + 2][i] = temp3;
                B[j + 3][i] = temp4;

                B[j][i + 4] = temp5;  //将A中右上 存到B中的右上，这也是全部命中的
                B[j + 1][i + 4] = temp6;
                B[j + 2][i + 4] = temp7;
                B[j + 3][i + 4] = temp8;
            }
            i -= 4;
            // 处理A中左下
            for (cnt = 0; cnt < 4; ++j, ++cnt) {  // 枚举0~8中的每一行，一行八列
                temp1 = A[i + 4][j];
                temp2 = A[i + 5][j];
                temp3 = A[i + 6][j];
                temp4 = A[i + 7][j];  // 拿到左下的元素

                // 因为我们这里本来就要处理 右上，所以这不会带来更多的miss
                temp5 = B[j][i + 4];  // 拿到我们之前赋给B右上的，也就是A右上的元素
                temp6 = B[j][i + 5];
                temp7 = B[j][i + 6];
                temp8 = B[j][i + 7];

                B[j][i + 4] = temp1;  //将左下翻转到B右上
                B[j][i + 5] = temp2;
                B[j][i + 6] = temp3;
                B[j][i + 7] = temp4;

                // 这一步，也不会带来更多的miss因为 i + 4 <= j + 4 <= i + 7 恒成立
                // 所以每次带来的 eviction 和 store B[j][i +
                // 4]只会带来一次MISS，和原来的操作是一样的

                B[j + 4][i] = temp5;  //将原B右上 赋值到 B左下， 这样A右上也就完成了翻转
                B[j + 4][i + 1] = temp6;
                B[j + 4][i + 2] = temp7;
                B[j + 4][i + 3] = temp8;
            }
            j -= 4;
            j += 4;  // 处理第四块 右下
            for (i += 4, cnt = 0; cnt < 4; ++cnt, ++i) {
                temp1 = A[i][j];  // 第四块没有任何改动， 和原来效果是一样的
                temp2 = A[i][j + 1];
                temp3 = A[i][j + 2];
                temp4 = A[i][j + 3];

                B[j][i] = temp1;
                B[j + 1][i] = temp2;
                B[j + 2][i] = temp3;
                B[j + 3][i] = temp4;
            }
            i -= 8, j -= 4;
        }
    }
  } else {
    for (i = 0; i < N; i += 16) {
      for (j = 0; j < M; j += 16) {
        for (ii = i; ii < i + 16 && ii < N; ii++) {
          for (jj = j; jj < j + 16 && jj < M; jj++) {
            a0 = A[ii][jj];
            B[jj][ii] = a0;
          }
        }
      }
    }
  }
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    // registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

