#include <stdio.h>
#include <math.h>
#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
int iabs(int a);
void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color);
void meshgrid (tgaImage *image, Model *model);
int sign (int a);
void triangle_2 (int x1, int x2, int x3, 
                 int y1, int y2, int y3, 
                 int z1, int z2, int z3, 
                 tgaImage *image, tgaColor color);
void meshgrid_2(tgaImage *image, Model *model);
void normal(int x1, int x2, int x3, 
            int y1, int y2, int y3, 
            int z1, int z2, int z3, Vec3* n);
double intensity(Vec3 light, Vec3 n);
void triangle_3 (int x1, int x2, int x3, 
                 int y1, int y2, int y3, 
                 int z1, int z2, int z3, 
                 tgaImage *image, tgaColor color, int (*zbuffer)[700][700]);
void meshgrid_3(tgaImage *image, Model *model); 

int main(int argc, char **argv){
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    
    Model *model = loadFromObj(argv[2]);
    int height = 700; 
    int width = 700;
    
    tgaImage * image = tgaNewImage(height, width, RGB);
    
    meshgrid_2(image, model);
   // meshgrid_3(image, model);

    if (-1 == tgaSaveToFile(image, argv[1])) {
        perror("tgaSateToFile");
        rv = -1;
    } 

    tgaFreeImage(image);
    freeModel(model);    
    return rv;
}

void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color){
    int steep = 0;
    if (iabs(y1 - y0) > iabs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }
    int dx = x1 - x0; 
    int dy = y1 - y0;
    int de = 2*iabs(dy);
    int e = 0;
    int y = y0;
    for (int x = x0; x <= x1; ++x) {
        if (steep == 1) {
            tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color); 
        } else {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
        }
        e = e + de;
        if (e > dx) {
            y = y + sign(dy); 
            e = e - 2*dx;
        }

    }
}

int sign (int a){
   int b;
   if (a > 0) {
      b = 1;
   }
   if (a == 0) {
      b = 0;
   }
   if (a < 0) {
      b = -1;
   }
   return b;
}


void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int iabs(int a) {
    return (a >= 0) ? a : -a;
}


void meshgrid(tgaImage *image, Model *model) {
    int i, j; 
    tgaColor white = tgaRGB(255, 255, 255);
    for (i = 0; i < model->nface; ++i) {
        int screen_coords[3][2]; // Переводим в экранные координаты
        for (j = 0; j < 3; ++j) {
            Vec3 *v = &(model->vertices[model->faces[i][3*j]]);
            screen_coords[j][0] = ((*v)[0] + 1) * image->width / 2;
            screen_coords[j][1] = (1 - (*v)[1]) * image->height / 2;
        }
        // Рисуем 3 ребра
        for (j = 0; j < 3; ++j) {
            line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j+1)%3][0], screen_coords[(j+1)%3][1], white);
        }
    }  
 }

void triangle_2 (int x1, int x2, int x3, 
                int y1, int y2, int y3, 
                int z1, int z2, int z3, 
                tgaImage *image, tgaColor color) {
     if (y1 == x2 && y1 == x3) return;
     if (y1 > y2) {
         swap(&y1, &y2);
         swap(&x1, &x2);
     }
     if (y1 > y3) {
         swap(&y1, &y3);
         swap(&x1, &x3);
     }
     if (y2 > y3) {
         swap(&y2, &y3);
         swap(&x2, &x3);
     }      
    int total_height = y3 - y1;
    for (int y = 0; y < total_height; y++) {
        _Bool second_half = y > y2 - y1 || y2 == y1;
        int segment_height = second_half ? y3 - y2 : y2 - y1;
        float alpha = (float)y/total_height;
        float beta  = (float)(y-(second_half ? y2 - y1 : 0))/segment_height;
        int xa = x1 + (x3 - x1)*alpha;
        int xb = second_half ? x2 + (x3 - x2)*beta : x1 + (x2 - x1)*beta;
        if (xa > xb){
            swap(&xa, &xb);
        }
        for (int x = xa; x <= xb; x++) {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)y1 + y, color);
        }
    }
}       

void meshgrid_2(tgaImage *image, Model *model) {
	int i, j;
    Vec3 light = {0, 0, -1};
	for (i = 0; i < model->nface; ++i) {
        double coords[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);            
            coords[j][0] = (*v)[0];
	  		coords[j][1] = (*v)[1];
            coords[j][2] = (*v)[2];       
		}        
        Vec3 n;
        normal(coords[0][0], coords[1][0], coords[2][0],
               coords[0][1], coords[1][1], coords[2][1],
               coords[0][2], coords[1][2], coords[2][2], &n);
        double I = intensity(light, n);
        if (I < 0){
             I = (-1) * I; 
             int sc[3][3];
            for (j = 0; j < 3; ++j){
		        sc[j][0] = (coords[j][0] + 1) * image->width / 2;
	            sc[j][1] = (1 - coords[j][1]) * image->height / 2;
                sc[j][2] = coords[j][2];    
            } 
            tgaColor color = tgaRGB(I*255, I*255, I*255);
            triangle_2(sc[0][0], sc[1][0], sc[2][0], 
                      sc[0][1], sc[1][1], sc[2][1], 
                      sc[0][2], sc[1][2], sc[2][2], image, color); 
        }                  
    }
}

void normal(int x1, int x2, int x3, 
            int y1, int y2, int y3, 
            int z1, int z2, int z3, Vec3* n){
    double a[3], b[3];
    a[0] = x2 - x1;
    a[1] = y2 - y1;
    a[2] = z2 - z1;
    b[0] = x3 - x1;
    b[1] = y3 - y1;
    b[2] = z3 - z1;
    (*n)[0] = a[1] * b[2] - b[1] * a[2];
    (*n)[1] = - a[0] * b[2] + b[0] * a[2];
    (*n)[2] = a[0] * b[1] - b[0] * a[1];
    double absn = sqrt((*n)[0] * (*n)[0] + (*n)[1] * (*n)[1] + (*n)[2] * (*n)[2]);
    (*n)[0] = (*n)[0] / absn;
    (*n)[1] = (*n)[1] / absn;
    (*n)[2] = (*n)[2] / absn;
}

double intensity(Vec3 light, Vec3 n){
    double I = 0.0;
    for (int i = 0; i < 3; ++i){
        I = I + light[i] * n[i];
    }
    return I;
}

void triangle_3 (int x1, int x2, int x3, 
                 int y1, int y2, int y3, 
                 int z1, int z2, int z3, 
                 tgaImage *image, tgaColor color, int (*zbuffer)[700][700]) {
     if (y1 == x2 && y1 == x3) return;
     if (y1 > y2) {
         swap(&y1, &y2);
         swap(&x1, &x2);
         swap(&z1, &z2);
     }
     if (y1 > y3) {
         swap(&y1, &y3);
         swap(&x1, &x3);
         swap(&z1, &z3);
     }
     if (y2 > y3) {
         swap(&y2, &y3);
         swap(&x2, &x3);
         swap(&z2, &z3);
     }      
    int total_height = y3 - y1;
    for (int y = 0; y < total_height; y++) {
        _Bool second_half = y > y2 - y1 || y2 == y1;
        int segment_height = second_half ? y3 - y2 : y2 - y1;
        float alpha = (float)y/total_height;
        float beta  = (float)(y-(second_half ? y2 - y1 : 0))/segment_height;
        int xa = x1 + (x3 - x1)*alpha;
        int za = z1 + (z3 - z1)*alpha;
        int xb = second_half ? x2 + (x3 - x2)*beta : x1 + (x2 - x1)*beta;
        int zb = second_half ? z2 + (z3 - z2)*beta : z1 + (z2 - z1)*beta;        
        if (xa > xb){
            swap(&xa, &xb);
            swap(&za, &zb);
        }
        for (int x = xa; x <= xb; x++) {
            float gamma = xa == xb ? 1.0 : (float)(x - xa)/(float)(xb - xa);
            int z = (float)za + (float)(zb - za)*gamma;
            if (z > (*zbuffer)[x][y1 + y]){
                  (*zbuffer)[x][y1 + y] = z;
                  tgaSetPixel(image, (unsigned int)x, (unsigned int)y1 + y, color);
            }
        }
    }
} 

void meshgrid_3(tgaImage *image, Model *model) {
	int i, j;
    Vec3 light = {0, 0, -1};
    int zbuffer[700][700];
    for (i = 0; i < 700; ++i){
         for (j = 0; j < 700; ++j){
             zbuffer[i][j] = -1;
          }
    }
	for (i = 0; i < model->nface; ++i) {
        double coords[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);            
            coords[j][0] = (*v)[0];
	  		coords[j][1] = (*v)[1];
            coords[j][2] = (*v)[2];       
		}        
        Vec3 n;
        normal(coords[0][0], coords[1][0], coords[2][0],
               coords[0][1], coords[1][1], coords[2][1],
               coords[0][2], coords[1][2], coords[2][2], &n);
        double I = intensity(light, n);
        if (I > 0) continue;
        I = -I; 
        int sc[3][3];
        for (j = 0; j < 3; ++j){
		    sc[j][0] = (coords[j][0] + 1) * image->width / 2;
	        sc[j][1] = (1 - coords[j][1]) * image->height / 2;
            sc[j][2] = coords[j][2];    
        } 
        tgaColor color = tgaRGB(I*255, I*255, I*255);
        triangle_3(sc[0][0], sc[1][0], sc[2][0], 
                   sc[0][1], sc[1][1], sc[2][1], 
                   sc[0][2], sc[1][2], sc[2][2], image, color, &zbuffer);         
    }
}
