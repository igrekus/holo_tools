//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dir.h>
#include <fftw3.h>
#include <conio.h>

#pragma hdrstop

#include "fnc.cpp"

#if !defined(LEPTONICA_ALLHEADERS_H)
#   include "lepton\allheaders.h"
#endif
//---------------------------------------------------------------------------

#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif

//#define KINO_SIZE_X 2048
//#define KINO_SIZE_Y 2048
//#define KINO_SIZE_X 1024
//#define KINO_SIZE_Y 1024
#define KINO_SIZE_X 512
#define KINO_SIZE_Y 512
//#define KINO_SIZE_X 50
//#define KINO_SIZE_Y 50

#define PLANE_RE 0
#define PLANE_IM 1

#define RND_TRUE 1
#define RND_FALSE 0

#define BIN_TRUE 1
#define BIN_FALSE 0

#pragma argsused

int _tmain(int argc, _TCHAR* argv[])
{
//   static fftw_complex com_image[KINO_SIZE_X][KINO_SIZE_Y];
//   static double double_amp[KINO_SIZE_X][KINO_SIZE_Y];

   double *double_amp;
   double *double_phs;
   unsigned char *bitmap_re;
//   unsigned char *bitmap_im;
   fftw_complex *com_image;

   fftw_plan plan_fwd, plan_bwd;
   int Nx=KINO_SIZE_X;
   int Ny=KINO_SIZE_Y;
   int bin=0;
   l_int32 i, j, w, h, bpp, format;
   l_uint32 val32;
//   l_uint32 *data;

   PIX *pix_in;
//   PIXCMAP *pix_in_cmap;

   char input_file_name[50];                           // Имя входного файла с изображением (.arx)
   char *output_file_name;                             // Имя выходного файла (.res)
   char *tmpstr;
   char fmtstr[20];

//   char output_im_file_name[20]="outim.res";      // Имя выходного файла (.res)
//   FILE *in, *out_re;                               // Идентификаторы входного и выходного файлов
//   FILE *out_im;

   if (argc == 1) {
      printf("\nUsage: kinview.exe <kinoform image file>\n\n8 BPP grayscale\\indexed color BMP, GIF, PNG, TIFF are supported.\n");
      //getch();
      return 1;
   } else {
      strcpy(input_file_name, argv[1]);
   }

//------------------ Считываем файл leptonica-ой ----------------------------
   printf("\nReading input...");
   pix_in=pixRead(input_file_name);

   if (!pix_in) {
      printf("Could not read input file.");
      return 1;
   }

   w=pix_in->w;
   h=pix_in->h;
   Nx=w;
   Ny=h;
   bpp=pix_in->d;
   format=pix_in->informat;

   printf("ok.\n");

   tmpstr=strdup(input_file_name);
   output_file_name=strtok(tmpstr, ".");

   switch (format) {
      case IFF_BMP: {                         // =1
         strcat(output_file_name, "-v.bmp\0");
//         fmtstr=strdup("BMP\0");
         strcpy(fmtstr, "BMP\0");
         break;
      }
      case IFF_PNG: {                         // =3
         strcat(output_file_name, "-v.png\0");
//         fmtstr=strdup("PNG\0");
         strcpy(fmtstr, "PNG\0");
         break;
      }
      case IFF_TIFF: {                        // =4
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("uncompressed TIFF\0");
         strcpy(fmtstr, "uncompressed TIFF\0");
         break;
      }
      case IFF_TIFF_PACKBITS: {               // =5
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("Packbits TIFF\0");
         strcpy(fmtstr, "packbits TIFF\0");
         break;
      }
      case IFF_TIFF_RLE: {                     // =6
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("RLE TIFF\0");
         strcpy(fmtstr, "RLE TIFF\0");
         break;
      }
      case IFF_TIFF_G3: {                     // =7
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("G3 TIFF\0");
         strcpy(fmtstr, "G3 TIFF\0");
         break;
      }
      case IFF_TIFF_G4: {                     // =8
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("G4 TIFF\0");
         strcpy(fmtstr, "G4 TIFF\0");
         break;
      }
      case IFF_TIFF_LZW: {                     // =9
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("LZW TIFF\0");
         strcpy(fmtstr, "LZW TIFF\0");
         break;
      }
      case IFF_TIFF_ZIP: {                     // =10
         strcat(output_file_name, "-v.tif\0");
//         fmtstr=strdup("LZW TIFF\0");
         strcpy(fmtstr, "ZIP TIFF\0");
         break;
      }
      case IFF_GIF: {                         // =13
         strcat(output_file_name, "-v.gif\0");
//         fmtstr=strdup("GIF\0");
         strcpy(fmtstr, "GIF\0");
         break;
      }
      default: {
         printf("\nUnsupported format file.\n");
         pixDestroy(&pix_in);
//         getch();
         return 2;
      }
   }

   printf("\nImage info: \n  H: %d \n  W: %d \nBPP: %d \nFMT: %s\n\n", w, h, bpp, fmtstr);
   printf("Output file: %s\n\n", output_file_name);
   if (bpp>8) {
      printf("Unsupported BPP.\n");
      pixDestroy(&pix_in);
      return 3;
   }

//   output_file_name = (char*) malloc(50*sizeof(char));

//   pix_in_cmap=pixGetColormap(pix_in);
//   data = pixGetData(pix_in);

//------------------ Выдеяем память под массивы -----------------------------
   double_amp = (double*) fftw_malloc(w*h*sizeof(double));
   double_phs = (double*) fftw_malloc(w*h*sizeof(double));
   bitmap_re = (unsigned char*) fftw_malloc(w*h*sizeof(unsigned char));
//   bitmap_im = (unsigned char*) malloc(w*h*sizeof(unsigned char));
   com_image = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*w*h);

   if (double_amp==NULL) {
      printf("Not enough memory.");
      return 1;
   };
   if (double_phs==NULL) {
      printf("Not enough memory.");
      return 1;
   };
   if (bitmap_re==NULL) {
      printf("Not enough memory.");
      return 1;
   };
/*   if (bitmap_im==NULL) {
      printf("Not enough memory.");
      return 1;
   };*/
   if (com_image==NULL) {
      printf("Not enough memory.");
      return 1;
   };

   for (i = 0; i < h; i++) {
      for (j = 0; j < w; j++) {
         pixGetPixel(pix_in, j, i, &val32);
         *(bitmap_re+j*w+i)=(unsigned char)val32;
//         printf("%d ", *(bitmap_re+j*w+i));
      }
   }

   printf("Retrieveing image...");

//---------- Подготовка FFTW ------------------------------------------------
   plan_fwd = fftw_plan_dft_2d(Nx, Ny, com_image, com_image, FFTW_FORWARD, FFTW_ESTIMATE);
   plan_bwd = fftw_plan_dft_2d(Nx, Ny, com_image, com_image, FFTW_BACKWARD, FFTW_ESTIMATE);
//   plan_fwd = fftw_plan_dft_2d(Nx, Ny, &com_image[0][0], &com_image[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
//   plan_bwd = fftw_plan_dft_2d(Nx, Ny, &com_image[0][0], &com_image[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

//---------- Считываем растр фазового рельефа -------------------------------
//   fread(bitmap_re, KINO_SIZE_X*KINO_SIZE_Y*sizeof(unsigned char), 1, in);
//---------- Преобразовываем растр в массив double и нормируем на 1 ---------
   bitmap_to_double(double_phs, bitmap_re, w, h);
   norm_double(double_phs, w, h, 1);
//---------- Бинарный рельеф? -----------------------------------------------
   bin=is_bin(bitmap_re, w, h);
   if (bin) {norm_pi_minus_pi(double_phs, w, h);} // Нормируем phs*2pi-pi
   else {norm_2pi_minus_pi(double_phs, w, h);}    // Нормируем phs*pi-pi
//---------- Амплитуда=1 ----------------------------------------------------
   double_fill_plane(double_amp, w, h, 1, RND_FALSE);
//---------- com_image=1*exp(i*phs) -----------------------------------------
   make_complex((double *)com_image, double_amp, double_phs, w, h);
//---------- Фурье +1 -------------------------------------------------------
   fftw_execute(plan_fwd);
//---------- Получаем амплитуды ---------------------------------------------
   com_get_amp((double *)com_image, double_amp, w, h);
//---------- Меняем четверти местами ----------------------------------------
   double_shift(double_amp, w, h);
//---------- Натягиваем палитру 0-255 grayscale -----------------------------
   double_rescale(double_amp, w, h, 0, 255);
//---------- Конвертим в битмап ---------------------------------------------
   double_to_bitmap(double_amp, bitmap_re, w, h);
//---------- Пишем в файл ---------------------------------------------------
//   fwrite(bitmap_re, KINO_SIZE_X*KINO_SIZE_Y*sizeof(unsigned char), 1, out_re);
   printf("ok.\n");
   printf("Writing output...");
   for (i = 0; i < h; i++) {
      for (j = 0; j < w; j++) {
         pixSetPixel(pix_in, j, i, (l_uint32)*(bitmap_re+j*w+i));
//         printf("%d ", (l_uint32)*(bitmap_re+j*w+i));
      }
   }

//   pixSetColormap(pix_in, pix_in_cmap);
//   pixAddGrayColormap8(pix_in);
   pixWrite(output_file_name, pix_in, format);
   printf("ok.\n");
//   com_echo(com_image, KINO_SIZE_X, KINO_SIZE_Y);
//   double_echo(double_amp, KINO_SIZE_X, KINO_SIZE_Y);

l1:
//--------- Завершаем работу с FFTW -----------------------------------------
   fftw_destroy_plan(plan_fwd);
   fftw_destroy_plan(plan_bwd);

//--------- Завершаем работу с PIX ------------------------------------------
   pixDestroy(&pix_in);

//---------- Освобождаем память ---------------------------------------------
   fftw_free(double_amp);
   fftw_free(double_phs);
   fftw_free(bitmap_re);
//   free(bitmap_im);
//   free(output_file_name);
//   free(tmpstr);
   fftw_free(com_image);
//   printf("\nDone.\n");
//   getch();
   return 0;
}
