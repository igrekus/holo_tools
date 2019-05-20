#pragma hdrstop
#pragma argsused

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dir.h>
#include <fftw3.h>
#include <conio.h>
#include <alloc.h>

#include "fnc.cpp"

#if !defined(LEPTONICA_ALLHEADERS_H)
#   include "lepton\allheaders.h"
#endif

#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif

#define KINO_SIZE_X 512
#define KINO_SIZE_Y 512

#define PLANE_RE 0
#define PLANE_IM 1
#define PLANE_BOTH 2

#define RND_TRUE 1
#define RND_FALSE 0

#define ABS_TRUE 1
#define ABS_FALSE 0

#define COLOR_MIN 0
#define COLOR_MAX 255

#define DEFAULT_CORRECTION 0.5

#define STRING_LENGTH 100 // Костыль для обхода косячащих имен файлов -- пофиксить!
#define FORMAT_LENGTH 10

int _tmain(int argc, _TCHAR* argv[])
{
//   static fftw_complex com_image[KINO_SIZE_X][KINO_SIZE_Y];
//   static double double_amp[KINO_SIZE_X][KINO_SIZE_Y];

   double *target_amp;
   double *source_amp;
   double *tmp_amp, *tmp_phs;

   double avg1, avg2, rmse;
   unsigned char *bitmap_re, *bitmap_im;
   fftw_complex *com_A, *com_B, *com_C, *com_D;

   fftw_plan plan_fwd, plan_bwd;
   int Nx=KINO_SIZE_X;
   int Ny=KINO_SIZE_Y;
   l_int32 i, j, w, h, bpp, format;
   l_uint32 val32;


   char input_file_name[STRING_LENGTH]="in.arx";       // Имя входного файла с изображением
   char *output_re_file_name;                          // Имя выходного файла для Re-части
   char *output_im_file_name;                          // Имя выходного файла для Im-части
   char *output_bin_file_name;                         // Имя выходного файла для бинаризованной Im-части
   char *tmpstr;
   char *fmtstr;
//   char dummystring[50]="thisisadummystring\0";

//   char source_file_name[20]="source.arx";             // Имя выходного файла
//   FILE *in, *out_re, *out_im, *src;                   // Идентификаторы входного и выходного файлов

   PIX *pix_in;

   if (argc == 1) {
      printf("\nUsage: kinprep.exe <target image file>\n\n8 BPP grayscale\\indexed color BMP, GIF, PNG, TIFF are supported.\n");
      return 1;
   } else {
      strcpy(input_file_name, argv[1]);
   }

//------------------ Считываем файл leptonica-ой ----------------------------
   printf("\nReading input...");
   pix_in=pixRead(input_file_name);

   if (!pix_in) {
      printf("Could not read input file.");
      return 2;
   }

   w=pix_in->w;
   h=pix_in->h;
   Nx=w;
   Ny=h;
   bpp=pix_in->d;
   format=pix_in->informat;

   printf("OK\n");

// Выделяем память под строки для имен файлов

   if ((output_re_file_name = (char *) malloc(sizeof(char)*STRING_LENGTH)) == NULL) {
       printf("Not enough memory to allocate buffer.\n");
       return 5;
   }
   if ((output_im_file_name = (char *) malloc(sizeof(char)*STRING_LENGTH)) == NULL) {
       printf("Not enough memory to allocate buffer.\n");
       return 5;
   }
   if ((output_bin_file_name = (char *) malloc(sizeof(char)*STRING_LENGTH)) == NULL) {
       printf("Not enough memory to allocate buffer.\n");
       return 5;
   }
   if ((tmpstr = (char *) malloc(sizeof(char)*STRING_LENGTH)) == NULL) {
       printf("Not enough memory to allocate buffer.\n");
       return 5;
   }
   if ((fmtstr = (char *) malloc(sizeof(char)*FORMAT_LENGTH)) == NULL) {
       printf("Not enough memory to allocate buffer.\n");
       return 5;
   }
      inline
   tmpstr=strdup(strtok(input_file_name, "."));
   strcat(tmpstr, "\0");

   memcpy(output_re_file_name, tmpstr, strlen(tmpstr));
   memcpy(output_im_file_name, tmpstr, strlen(tmpstr));
   memcpy(output_bin_file_name, tmpstr, strlen(tmpstr));

   switch (format) {
      case IFF_BMP: {                         // =1
         strcat(output_im_file_name, "-i.bmp\0");
         strcat(output_re_file_name, "-r.bmp\0");
         strcat(output_bin_file_name, "-ib.bmp\0");
//         fmtstr=strdup("BMP\0");
         strcpy(fmtstr, "BMP\0");
         break;
      }
      case IFF_PNG: {                         // =3
         strcat(output_im_file_name, "-i.png\0");
         strcat(output_re_file_name, "-r.png\0");
         strcat(output_bin_file_name, "-ib.png\0");
//         fmtstr=strdup("PNG\0");
         strcpy(fmtstr, "PNG\0");
         break;
      }
      case IFF_TIFF:                          // =4
      case IFF_TIFF_PACKBITS:                 // =5
      case IFF_TIFF_RLE:                      // =6
      case IFF_TIFF_G3:                       // =7
      case IFF_TIFF_G4:                       // =8
      case IFF_TIFF_LZW:                      // =9
      case IFF_TIFF_Z	IP: {                    // =10
         strcat(output_im_file_name, "-i.tif\0");
         strcat(output_re_file_name, "-r.tif\0");
         strcat(output_bin_file_name, "-ib.tif\0");
//         fmtstr=strdup("TIFF\0");
         strcpy(fmtstr, "TIFF\0");
         break;
      }
      case IFF_GIF: {                         // =13
         strcat(output_im_file_name, "-i.gif\0");
         strcat(output_re_file_name, "-r.gif\0");
         strcat(output_bin_file_name, "-ib.gif\0");
//         fmtstr=strdup("uncompressed TIFF\0");
         strcpy(fmtstr, "GIF\0");
         break;
      }
      default: {
         printf("\nUnsupported format file.\n");
         pixDestroy(&pix_in);
//         getch();
         return 3;
      }
   }

   printf("\nImage info: \n  W: %d \n  H: %d \nBPP: %d \nFMT: %s\n\n", w, h, bpp, fmtstr);

   if (bpp>8) {
      printf("Unsupported BPP.\n");
      pixDestroy(&pix_in);
      return 4;
   }

//   printf("Output re: %s\nOutput im: %s\n\n", output_re_file_name, output_im_file_name);
   printf("Output gray: %s\nOutput bin: %s\n\n", output_im_file_name, output_bin_file_name);

//------------------ Выдеяем память под массивы -----------------------------
   target_amp = (double*) fftw_malloc(Nx*Ny*sizeof(double));
   source_amp = (double*) fftw_malloc(Nx*Ny*sizeof(double));
   tmp_amp = (double*) fftw_malloc(Nx*Ny*sizeof(double));
   tmp_phs = (double*) fftw_malloc(Nx*Ny*sizeof(double));

   bitmap_re = (unsigned char*) fftw_malloc(Nx*Ny*sizeof(unsigned char));
   bitmap_im = (unsigned char*) fftw_malloc(Nx*Ny*sizeof(unsigned char));
   com_A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
   com_B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
   com_C = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
   com_D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

   if (target_amp==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (source_amp==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (bitmap_re==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (bitmap_im==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (com_A==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (com_B==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (com_C==NULL) {
      printf("Not enough memory.");
      return 5;
   };
   if (com_D==NULL) {
      printf("Not enough memory.");
      return 5;
   };

//---------- Подготовка FFTW ------------------------------------------------
   plan_fwd = fftw_plan_dft_2d(Nx, Ny, com_B, com_C, FFTW_FORWARD, FFTW_ESTIMATE);
   plan_bwd = fftw_plan_dft_2d(Nx, Ny, com_D, com_A, FFTW_BACKWARD, FFTW_ESTIMATE);
//   plan_fwd = fftw_plan_dft_2d(Nx, Ny, &com_I[0][0], &com_I[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
//   plan_bwd = fftw_plan_dft_2d(Nx, Ny, &com_I[0][0], &com_I[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

/*
Calculate modulus of pixel (complex) value
     NumberType abs = sqrt(real*real+imag*imag);

energy conservation scale factor
     NumberType factor = sqrt(sum*1.0)/(size*1.0);

with conservation
     if(__is_zero(abs))
     {
         real = (real/abs)*factor;
         imag = (imag/abs)*factor;
     }
     else
     {
         real = 0;
         imag = 0;
     }
     d_devPtr[x].x=real;
     d_devPtr[x].y=imag;
 }
//-------
twilight error

     d_error[x].x = 0;
     d_error[x].y = 0;

computed amplitude for current pixel
     T a_ampl = d_devPtr[x];
     NumberType abs_a = (a_ampl.x*a_ampl.x + a_ampl.y*a_ampl.y);

original amplitude for current pixel
     T A_ampl = d_original[x];
     NumberType abs_A = (A_ampl.x*A_ampl.x + A_ampl.y*A_ampl.y);

error margins
     NumberType t_light = 0.1;
     NumberType t_dark = 3.333e-4;
//
     NumberType err_light = (abs(abs_A-abs_a)/abs_A - t_light);
     NumberType err_dark = (abs_a-t_dark);
     if(A_ampl.x == 0)
     {//error for dark pixel
         if(err_dark>0)
         {
             d_error[x].x = err_dark;
         }
     }
     else
     {//error for light pixel
         if(err_light>0)
         {
             d_error[x].x= err_light*(t_dark/t_light);
         }
     }
 }

*/
//---------- Подготавливаем source_amp --------------------------------------
//   fread(bitmap_re, Nx*Ny*sizeof(unsigned char), 1, src);
//   bitmap_to_double(source_amp, bitmap_re, Nx, Ny);
//   norm_double(source_amp, Nx, Ny, 1);

//////////// Алгоритм Gerschberg-Saxton //////////////////////////////////////

//---------- Считываем растр для расчёта ------------------------------------
//   fread(bitmap_re, Nx*Ny*sizeof(unsigned char), 1, in);

//---------- Заполняем комплексную матрицу изображения ----------------------
//---------- re-часть - растр, ----------------------------------------------
   for (i = 0; i < h; i++) {
      for (j = 0; j < w; j++) {
         pixGetPixel(pix_in, j, i, &val32);
         *(bitmap_re+j*w+i)=(unsigned char)val32;
//         printf("%d ", *(bitmap_re+j*w+i));
      }
   }

   printf("Interating...\n\n");

   bitmap_to_double(target_amp, bitmap_re, Nx, Ny);
//---------- Нормируем re-часть на 1 ----------------------------------------
   norm_double(target_amp, Nx, Ny, 1);
//---------- Корректируем изображение по sinc -------------------------------
   double_sinc(target_amp, Nx, Ny, DEFAULT_CORRECTION);

//   memcpy(tmp_phs, target_amp, Nx*Ny*sizeof(double));
//goto l2;

//---------- Меняем четверти матрицы местами по диагонали -------------------
   double_shift(target_amp, Nx, Ny);
   set_plane((double *)com_D, target_amp, Nx, Ny, PLANE_RE);
//---------- im-часть - случайный шум ---------------------------------------
   com_fill_plane((double *)com_D, Nx, Ny, PLANE_IM, 255, RND_TRUE);
//   com_fill_plane((double *)com_D, Nx, Ny, PLANE_IM, 0, RND_FALSE);

//---------- Вычисляем исходную среднюю интенсивность -----------------------
//   avg1=double_avg(target_amp, Nx, Ny, ABS_TRUE);

   fftw_execute(plan_bwd);           // A = IFT(D) (D = target)

   for (i=0; i<9; i++) {
      printf("Pass: %d...", i+1);
      com_get_phase((double *)com_A, tmp_phs, Nx, Ny);
      double_fill_plane(tmp_amp, Nx, Ny, 1, RND_FALSE);
      make_complex((double *)com_B, tmp_amp, tmp_phs, Nx, Ny);

      fftw_execute(plan_fwd);        // C = FT(B)
      norm_fft((double *)com_C, Nx, Ny);

      com_get_phase((double *)com_C, tmp_phs, Nx, Ny);
      make_complex((double *)com_D, target_amp, tmp_phs, Nx, Ny);

      fftw_execute(plan_bwd);        // A = IFT(D)
      printf("OK\n");
//      writefile((double *)com_D);
   }

   com_get_phase((double *)com_A, tmp_phs, Nx, Ny);
   double_rescale(tmp_phs, Nx, Ny, COLOR_MIN, COLOR_MAX);

//---------- Пишем серый киноформ в файл ---------------------------------------
   double_to_bitmap(tmp_phs, bitmap_im, Nx, Ny);

   printf("\nWriting gray...");
   for (i = 0; i < h; i++) {
      for (j = 0; j < w; j++) {
         pixSetPixel(pix_in, j, i, (l_uint32)*(bitmap_im+j*w+i));
//         printf("%d ", (l_uint32)*(bitmap_re+j*w+i));
      }
   }
   pixWrite(output_im_file_name, pix_in, format);
   printf("OK\n");

   printf("\nBinarizing...\n\n");

   for (i=0; i<9; i++) {
      printf("Pass: %d...", i+1);

      com_get_phase((double *)com_A, tmp_phs, Nx, Ny);

      double_binarize(tmp_phs, Nx, Ny);
//      double_echo(tmp_phs, Nx, Ny);

      double_fill_plane(tmp_amp, Nx, Ny, 1, RND_FALSE);
      make_complex((double *)com_B, tmp_amp, tmp_phs, Nx, Ny);

      fftw_execute(plan_fwd);        // C = FT(B)
      norm_fft((double *)com_C, Nx, Ny);

      com_get_phase((double *)com_C, tmp_phs, Nx, Ny);
      make_complex((double *)com_D, target_amp, tmp_phs, Nx, Ny);

      fftw_execute(plan_bwd);        // A = IFT(D)

      printf("OK\n");
   }
/**/
/*
// Среднеквадратичное отклонение между модулем I2 и re-частью I
      get_plane_abs((double *)com_I2, double_amp, Nx, Ny, PLANE_RE);
      get_plane((double *)com_I, double_phs, Nx, Ny, PLANE_RE);
      rmse=com_rmse(double_amp, double_phs, Nx, Ny);

      printf("Pass: %d, RMSE: %2.4f\n", i+1, rmse);
   }
*/

//---------- Бинарный киноформ -----------
   com_get_phase((double *)com_A, tmp_phs, Nx, Ny);
//   double_echo(tmp_phs, Nx, Ny);

   double_binarize(tmp_phs, Nx, Ny);
l2:
   double_rescale(tmp_phs, Nx, Ny, COLOR_MIN, COLOR_MAX);

   double_to_bitmap(tmp_phs, bitmap_im, Nx, Ny);

   printf("\nWriting bin...");
   for (i = 0; i < h; i++) {
      for (j = 0; j < w; j++) {
         pixSetPixel(pix_in, j, i, (l_uint32)*(bitmap_im+j*w+i));
//         printf("%d ", (l_uint32)*(bitmap_re+j*w+i));
      }
   }
   pixWrite(output_bin_file_name, pix_in, format);
   printf("OK\n");

l1:
//--------- Завершаем работу с FFTW -----------------------------------------
   fftw_destroy_plan(plan_fwd);
   fftw_destroy_plan(plan_bwd);
//--------- Завершаем работу с PIX ------------------------------------------
   pixDestroy(&pix_in);
//---------- Освобождаем память ---------------------------------------------

   free(output_re_file_name);
   free(output_im_file_name);
   free(output_bin_file_name);
   free(tmpstr);
   free(fmtstr);

   fftw_free(target_amp);
   fftw_free(source_amp);
   fftw_free(tmp_amp);
   fftw_free(tmp_phs);
   fftw_free(bitmap_re);
   fftw_free(bitmap_im);
//   free(tmpstr);
   fftw_free(com_A);
   fftw_free(com_B);
   fftw_free(com_C);
   fftw_free(com_D);

   printf("\nDone.\n");
//   getch();
   return 0;
}
//---------------------------------------------------------------------------
