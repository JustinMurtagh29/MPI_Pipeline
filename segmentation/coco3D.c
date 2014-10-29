/**
 * \file coco3D.c
 * \brief Fast Connected Components in 3D.
 * This algorithm finds the connected components in 3D affinity maps.
 * @param cnn_x, cnn_y, cnn_y
 * @param neighbourhood
 * @return segmented-image
 */

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

int NX, NY, NZ;
int labels[16384][2]; // store equiv. class and rank
int max_label = -1;

/*
void consolidate()
{
  int i, k;
  int *new_labels = (int *) mxMalloc(max_labels * sizeof(int));
  int nlc = 0; // new label counter
  for (i = 0; i < max_label; i++) {
    if (labels[i][0] < 0) continue;
    new_labels[nlc] = labels[i][0];
    for (k = i; k < max_label; k++) {
      if (labels[k][0] == new_labels[nlc]) {
	labels[k][0] = -nlc;
      }
    }
    nlc++;
  }
  for (i = 0; i < max_label; i++) {
    labels[i][0] = -labels[i][0];
  }
  
}
*/

void consolidate(int *result)
{
  int i, k;
  int nlc = -1; // new label counter
  for (i = 0; i < NX*NY*NZ; i++) {
    if (result[i] <= 0) continue;
    for (k = i+1; k < NX*NY*NZ; k++) {
      if (result[k] == result[i]) {
	result[k] = nlc;
      }
    }
    result[i] = nlc;
    nlc--;
  }
  for (i = 0; i < NX*NY*NZ; i++) {
    result[i] = -result[i];
  }
  
}

unsigned char getmin(unsigned char a, unsigned char b)
{
    if (a == 0 && b == 0) return 255;
    if (a == 0) return b;
    if (b == 0) return a;
    if (a > b) return b; else return a;
}

int find(int x)
{
    while (x != labels[x][0]) {
        labels[x][0] = labels[labels[x][0]][0];
        labels[x][1] = labels[labels[x][0]][1];
        x = labels[x][0];
    }
    return x;
}

void unite(int x, int y)
{
    int xRoot, yRoot;
    xRoot = find(x);
    yRoot = find(y);
    if (xRoot == yRoot) return;
    if (labels[xRoot][1] < labels[yRoot][1]) {
        labels[xRoot][0] = yRoot;
    } else if (labels[xRoot][1] > labels[yRoot][1]) {
        labels[yRoot][0] = xRoot;
    } else {
        labels[yRoot][0] = xRoot;
        labels[xRoot][1] = labels[xRoot][1] + 1;
    }
}

void coco3D(int *result, unsigned char *cx, unsigned char *cy, unsigned char *cz, int nb_mode)
{
    int i;
    int x, y, z;
    unsigned char c;
    unsigned char *nb;
    unsigned char mini;
    unsigned char data;
    int lab = 1;

    if (nb_mode != 26 && nb_mode != 18 && nb_mode != 6) {
        nb_mode = 6;
    }
    // fix here
    nb = (unsigned char *) mxMalloc(nb_mode * sizeof(unsigned char));

    for (x = 0; x < 256; x++) {
        labels[x][0] = 0;
        labels[x][1] = 0;
    }

    for (z = 0; z < NZ; z++) {
        for (y = 0; y < NY; y++) {
            for (x = 0; x < NX; x++) {
                data = (cx[z*NX*NY + y*NX + x] + cy[z*NX*NY + y*NX + x] + cz[z*NX*NY + y*NX + x]) / 2;
                if (data == 0) {
                    cx[z*NX*NY + y*NX + x] = 0;
                    cy[z*NX*NY + y*NX + x] = 0;
                    cz[z*NX*NY + y*NX + x] = 0;
                    continue;
                }

                if (x == 0) nb[0] = 0; else {
                    nb[0] = cx[z*NX*NY + y*NX + x-1];
                }
                if (y == 0) nb[1] = 0; else {
                    nb[1] = cy[z*NX*NY + (y-1)*NX + x];
                }
                if (z == 0) nb[2] = 0; else {
                    nb[2] = cz[(z-1)*NX*NY + y*NX + x];
                }

                if (nb[0] == 0 && nb[1] == 0 && nb[2] == 0) {
                    // if all are 0, then create a new label
                    c = lab;
                    labels[lab][0] = lab;
                    labels[lab][1] = 0;
                    lab++;
                } else {
                    // choose the min label and store the equivalence relationship
                    mini = getmin(getmin(nb[0], nb[1]), nb[2]);
                    for (i = 0; i < 3; i++) {
                        if (nb[i] != 0 && nb[i] != mini) {
                            unite(nb[i], mini);
                        }
                    }
                    c = mini;
                }
                cx[z*NX*NY + y*NX + x] = c;
                cy[z*NX*NY + y*NX + x] = c;
                cz[z*NX*NY + y*NX + x] = c;
            }
        }
    }

    max_label = lab;
    printf("max label: %d\n", max_label);
    printf("%d %d %d\n", NX, NY, NZ);

    printf("\n");
    printf("label: ");
    for (x = 0; x < lab; x++) {
        printf("%d ", x);
    }
    printf("\n\n");
    
    printf("equiv: ");
    for (x = 0; x < lab; x++) {
        printf("%d ", labels[x][0]);
    }
    printf("\n\n");

    // second pass
    for (z = 0; z < NZ; z++) {
        for (y = 0; y < NY; y++) {
            for (x = 0; x < NX; x++) {
	        // cx is taken, because now cx = cy = cz
                result[z*NX*NY + y*NX + x] = find(cx[z*NX*NY + y*NX + x]);
            }
        }
    }

    consolidate(result);
}


/**
 * mexFunction
 * Matlab routine should be called like this:
 * result = coco3D(cnn_x, cnn_y, cnn_z);
 * cnn_x, cnn_y and cnn_z are images created by the CNN
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int x, y, z;
    unsigned char *cx;
    unsigned char *cy;
    unsigned char *cz;
    int *result;
    const mwSize *dim;

    dim = mxGetDimensions(prhs[0]);
    NX = dim[0];
    NY = dim[1];
    NZ = dim[2];
    printf("%d %d %d\n", NX, NY, NZ);

    cx = (unsigned char *) mxGetData(prhs[0]);
    cy = (unsigned char *) mxGetData(prhs[1]);
    cz = (unsigned char *) mxGetData(prhs[2]);

    plhs[0] = (mxArray *) mxCreateNumericArray(3, dim, mxINT32_CLASS, mxREAL);
    result = (int *) mxGetData(plhs[0]);

    coco3D(result, cx, cy, cz, 6);

    printf("DONE\n");
}

