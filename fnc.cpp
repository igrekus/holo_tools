//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dir.h>
#include <fftw3.h>
#include <time.h>
#include <math.h>

#pragma hdrstop

#include "fnc.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#define PLANE_RE 0
#define PLANE_IM 1
#define PLANE_BOTH 2

#define RND_TRUE 1
#define RND_FALSE 0

#define ABS_TRUE 1
#define ABS_FALSE 0

#define COLOR_MIN 0
#define COLOR_MAX 255

int bitmap_to_double(double *double_buff, unsigned char *bmp_buff, int size_x, int size_y) {
// Конвертим растр в массив double
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=(double)(*(bmp_buff+j*size_x+i));
//         printf(" %2.2f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int com_bitmap_to_plane(double *com_buff, unsigned char *bmp_buff, int size_x, int size_y, int plane) {
// Заполняем одну из плоскостей комплексной матрицей заданным растром
// В программе расшифровки изображенией эта функция работает по-другому!
// plane - действительная/мнимая (0/1) части комплексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(com_buff+2*(j*size_x+i)+plane)=(double)(*(bmp_buff+j*size_x+i));
//         printf(" %2.2f", *(com_buff+2*(j*size_x+i)+plane));
      }
   }
   return 0;
}

int com_copy_image(double *com_src, double *com_dst, int size_x, int size_y) {
// Копировать src комплексную матрицу в dst матрицу
// todo Оптимизировать под линейный массив размером size_x*size_y
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(com_dst+2*(j*size_x+i)+PLANE_RE)=*(com_src+2*(j*size_x+i)+PLANE_RE);
         *(com_dst+2*(j*size_x+i)+PLANE_IM)=*(com_src+2*(j*size_x+i)+PLANE_IM);
//         printf(" %2.2f+%2.2f", *(com_dst+2*(j*size_x+i)+PLANE_RE), *(com_dst+2*(j*size_x+i)+PLANE_IM));
      }
   }
   return 0;
}

int com_echo(double *com_buff, int size_x, int size_y) {
// Вывод содержимого комплексного массива на экран для отладки
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         printf("%2.4f + %2.4f ", *(com_buff+2*(j*size_x+i)+PLANE_RE), *(com_buff+2*(j*size_x+i)+PLANE_IM));
      }
   }
   return 0;
}

int com_fill_plane(double *com_buff, int size_x, int size_y, int plane, int  val, int rnd) {
// Заполняем заданную плоскость заданным значением
// plane - действительная/мнимая = (0)/(1) части матрицы
   int i, j;
   if (rnd) {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(com_buff+2*(j*size_x+i)+plane)=random(val);
//            printf(" %2.2f", *(com_buff+2*(j*size_x+i)+plane));
         }
      }
   } else {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(com_buff+2*(j*size_x+i)+plane)=val;
//            printf(" %2.2f", *(com_buff+2*(j*size_x+i)+plane));
         }
      }
   }
   return 0;
}

int com_get_amp(double *com_buff, double *double_buff, int size_x, int size_y) {
// Получаем амплитуды копмлексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=sqrt(pow(*(com_buff+2*(j*size_x+i)+PLANE_RE),2)+pow(*(com_buff+2*(j*size_x+i)+PLANE_IM),2));
//         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int com_get_phase(double *com_buff, double *double_buff, int size_x, int size_y) {
// Получаем амплитуды копмлексной матрицы
// Не совпадают минусы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=atan2(*(com_buff+2*(j*size_x+i)+PLANE_IM), *(com_buff+2*(j*size_x+i)+PLANE_RE));
//         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int com_plane_to_bitmap(double *com_buff, unsigned char *bmp_buff, int size_x, int size_y, int plane) {
// Выдираем одну из частей комплексной матрицы и конвертим её в битмап
// plane - действительная/мнимая (0/1) части комплексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(bmp_buff+j*size_x+i)=*(com_buff+2*(j*size_x+i)+plane);
      }
   }
   return 0;
}

double com_rmse(double *double_buff1, double *double_buff2, int size_x, int size_y) {
// Вычисляем среднеквадратичное отклонение между buff1 и buff2
// Квадратный корень из среднего по всем (buff1-buff2)^2
   int i, j;
   double tmp=0;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         tmp=tmp+pow(*(double_buff1+j*size_x+i)-*(double_buff2+j*size_x+i),2);
      }
   }
   tmp=sqrt(tmp/(size_x*size_y));
   return tmp;
}

int com_shift(double *com_buff, int size_x, int size_y) {
   div_t x, y;
   int i, j;
   int nx, ny;
   double tmp_re, tmp_im;
   x=div(size_x,2);
   y=div(size_y,2);
   nx=size_x/2;
   ny=size_y/2;
   if ((x.rem==0)&&(y.rem==0)) {
      for (j=0; j<ny; j++) {
         for (i=0; i<nx; i++) {
            tmp_re=*(com_buff+2*(j*size_x+i)+PLANE_RE);
            tmp_im=*(com_buff+2*(j*size_x+i)+PLANE_IM);
            *(com_buff+2*(j*size_x+i)+PLANE_RE)=*(com_buff+2*((ny+j)*size_x+(nx+i))+PLANE_RE);
            *(com_buff+2*(j*size_x+i)+PLANE_IM)=*(com_buff+2*((ny+j)*size_x+(nx+i))+PLANE_IM);
            *(com_buff+2*((ny+j)*size_x+(nx+i))+PLANE_RE)=tmp_re;
            *(com_buff+2*((ny+j)*size_x+(nx+i))+PLANE_IM)=tmp_im;

            tmp_re=*(com_buff+2*(j*size_x+(nx+i))+PLANE_RE);
            tmp_im=*(com_buff+2*(j*size_x+(nx+i))+PLANE_IM);
            *(com_buff+2*(j*size_x+(nx+i))+PLANE_RE)=*(com_buff+2*((ny+j)*size_x+i)+PLANE_RE);
            *(com_buff+2*(j*size_x+(nx+i))+PLANE_IM)=*(com_buff+2*((ny+j)*size_x+i)+PLANE_IM);
            *(com_buff+2*((ny+j)*size_x+i)+PLANE_RE)=tmp_re;
            *(com_buff+2*((ny+j)*size_x+i)+PLANE_IM)=tmp_im;
         }
      }
   } else {    // пока не работает
   }
   return 0;
}

double double_avg(double *double_buff, int size_x, int size_y, int absol) {
// Вычисляем среднее по вещественной матрице
   int i, j;
   double tmp=0;
   if (absol) {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            tmp=tmp+fabs(*(double_buff+j*size_x+i));
//            printf(" %2.4f", *(double_buff+j*size_x+i));
         }
      }
   } else {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            tmp=tmp+*(double_buff+j*size_x+i);
         }
      }
   }
   tmp=tmp/(size_x*size_y);
   return tmp;
}

// Бинаризация по углам - приемлемо
int double_binarize(double *double_buff, int size_x, int size_y) {
// Бинаризуем массив double
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         if ((*(double_buff+j*size_x+i))>0) *(double_buff+j*size_x+i)=M_PI;
         if ((*(double_buff+j*size_x+i))<=0) *(double_buff+j*size_x+i)=0;
      }
   }
   return 0;
}
/**/

/*// Бинаризация по среднему уровню - приемлемо
int double_binarize(double *double_buff, int size_x, int size_y) {
// Бинаризуем массив double
   int i, j;
   double mid;
   double min=*(double_buff+0*size_x+0);
   double max=*(double_buff+0*size_x+0);

   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         if (min>(*(double_buff+j*size_x+i))) min=*(double_buff+j*size_x+i);
         if (max<(*(double_buff+j*size_x+i))) max=*(double_buff+j*size_x+i);
      }
   }

   mid=(max-min)/2;

   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         if ((*(double_buff+j*size_x+i))>=mid) {
            *(double_buff+j*size_x+i)=255;
         } else {
            *(double_buff+j*size_x+i)=0;
         }
      }
   }

   return 0;
}
/**/

/*// Бинаризация методом Otsu - плохой результат
int double_binarize(double *double_buff, int size_x, int size_y) {
// Бинаризуем массив double
   int i, j;
   double min=*(double_buff+0*size_x+0);
   double max=*(double_buff+0*size_x+0);
   double temp, temp1;
   int *hist;
   int hist_size;
   int temp_ind=0;
   double alpha, beta, threshold=0;
   double sigma, max_sigma=-1;
   double w1, a;

// Выбираем порог бинаризации методом Otsu
   double_rescale(double_buff, size_x, size_y, COLOR_MIN, COLOR_MAX);

   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         temp=*(double_buff+j*size_x+i);
         if (temp<min) {min=temp;}
         if (temp>max) {max=temp;}
      }
   }

// Память под гистограмму
   hist_size=max-min+1;

   if((hist=(int*)malloc(sizeof(int)*hist_size))==NULL) return -1;

   for (i=0; i<hist_size; i++) {
      hist[i]=0;
   }

// Строим гистограмму
   for (i=0; i<size_x*size_y; i++) {
      temp_ind=(int)(*(double_buff+i))-min;
      hist[temp_ind]++;
   }

   temp=temp1=0;
   alpha=beta=0;

// Для расяета матожидания 1го класса

   for (i=0; i<=(max-min); i++) {
      temp=temp+i*hist[i];
      temp1=temp1+hist[i];
   }

// Основной цикл поиска порога
// Проверяем все полутона, ищем то, при котором внутриклассовая дисперсия -> min

   for (i=0; i<=(max-min); i++) {
      alpha=alpha+i*hist[i];
      beta=beta+hist[i];

      w1=(double)beta/temp1;
      a=(double)alpha/beta-(double)(temp-alpha)/(temp-alpha);
      sigma=w1*(1-w1)*a*a;

      if (sigma>max_sigma) {
         max_sigma=sigma;
         threshold=i;
      }
   }
   free(hist);

   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         if ((*(double_buff+j*size_x+i))>=threshold) {
//            *(double_buff+j*size_x+i)=255;
            *(double_buff+j*size_x+i)=255;
         } else {
            *(double_buff+j*size_x+i)=0;
         }
//         printf("%2.2f ", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}
/**/

int double_echo(double *double_buff, int size_x, int size_y) {
// Вывод содержимого массива doble на экран для отладки
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int double_fill_plane(double *double_buff, int size_x, int size_y, int val, int rnd) {
// Заполняем массив double заданным значением
// plane - действительная/мнимая = (0)/(1) части матрицы
   int i, j;
   if (rnd) {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(double_buff+j*size_x+i)=random(val);
         }
      }
   } else {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(double_buff+j*size_x+i)=val;
//            printf(" %2.2f", *(double_buff+j*size_x+i));
         }
      }
   }
   return 0;
}

int double_rescale(double *double_buff, int size_x, int size_y, int c_min, int c_max) {
// Здесь нормализуем массив на отрезке 0..255 для вывода
// grad=1 - плавный градиент от 0 да 255
// grad=0 - отрезаем все, что <0 (бинаризация)
   int i, j;
   double low=255;
   double high=-255;
   double period;
   double mid;

   for (j=0; j<size_y; j++) {
		for (i=0; i<size_x; i++) {
      	if ((*(double_buff+j*size_x+i))<low) {
            low=*(double_buff+j*size_x+i);
         }
      	if ((*(double_buff+j*size_x+i))>high) {
            high=*(double_buff+j*size_x+i);
         }
      }
   }
   period=(high-low)/(c_max-c_min);
   for (j=0; j<size_y; j++) {
    	for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=floor(((*(double_buff+j*size_x+i))-low)/period);
      }
   }

/* бинаризация
      mid=low+(high-low)/2;
      for (j=0; j<size_y; j++) {
	   	for (i=0; i<size_x; i++) {
            if (*(double_buff+j*size_x+i)<=mid_re) {
               *(double_buff+j*size_x+i)=0;
            } else {
               *(double_buff+j*size_x+i)=255;
            }
         }
      }
   }
*/
   return 0;
}

int double_shift(double *double_buff, int size_x, int size_y) {
// Меняем четверти массива местами по диагонали
   div_t x, y;
   int i, j;
   int nx, ny;
   double tmp;
   x=div(size_x,2);
   y=div(size_y,2);
   nx=size_x/2;
   ny=size_y/2;
   if ((x.rem==0)&&(y.rem==0)) {
      for (j=0; j<ny; j++) {
         for (i=0; i<nx; i++) {
            tmp=*(double_buff+j*size_x+i);
            *(double_buff+j*size_x+i)=*(double_buff+(ny+j)*size_x+(nx+i));
            *(double_buff+(ny+j)*size_x+(nx+i))=tmp;
            tmp=*(double_buff+j*size_x+(nx+i));
            *(double_buff+j*size_x+(nx+i))=*(double_buff+(ny+j)*size_x+i);
            *(double_buff+(ny+j)*size_x+i)=tmp;
         }
      }
   } else {    // пока не работает
      if ((x.rem==1)&&(y.rem==1)) {
         for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
               tmp=*(double_buff+j*size_x+i);
               *(double_buff+j*size_x+i)=*(double_buff+(ny+1+j)*size_x+(nx+1+i));
               *(double_buff+(ny+1+j)*size_x+(nx+1+i))=tmp;
               tmp=*(double_buff+j*size_x+(nx+1+i));
               *(double_buff+j*size_x+(nx+1+i))=*(double_buff+(ny+1+j)*size_x+i);
               *(double_buff+(ny+1+j)*size_x+i)=tmp;
            }
         }
      }
   }
   return 0;
}

int double_sinc(double *double_buff, int size_x, int size_y, double corr) {
// Корректируем интенсивность входного изображения по sinc=sin(x)/x
   int i, j, N, r;
   double x, cs;
   double *sinc_array;

   N=size_y;
   cs=corr;

// Выделяем память под массив коэффициентов sinc
   sinc_array=(double*)fftw_malloc(1.5*N*sizeof(double));

   sinc_array[0]=1;

   for (i=1; i<N; i++) {
      x=i*M_PI/2/N;
      *(sinc_array+i)=cs*pow((sin(x)/x),2);
   }
//        sink2(1)=(2./PI)**2
   for (j=1; j<size_y; j++) {
      for (i=1; i<size_x; i++) {
         r=sqrt(pow(size_x/2-i,2)+pow(size_y/2-j,2));
//         printf(" r%d i%d j%d", r, i, j);
         *(double_buff+j*size_x+i)=*(double_buff+j*size_x+i)/sinc_array[r]/sinc_array[r];
      }
   }
   fftw_free(sinc_array);
   return 0;
}


int double_to_bitmap(double *double_buff, unsigned char *bmp_buff, int size_x, int size_y) {
// Конвертим массив double в битмап
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(bmp_buff+j*size_x+i)=*(double_buff+j*size_x+i);
      }
   }
   return 0;
}

int fill_double(double *double_buff, int size_x, int size_y, int val, int rnd) {
// Заполняем массив double заданным значением
   int i, j;
   if (rnd) {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(double_buff+j*size_x+i)=random(val);
//            printf(" %2.2f", *(double_buff+j*size_x+i));
         }
      }
   } else {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(double_buff+j*size_x+i)=val;
//            printf(" %2.2f", *(double_buff+j*size_x+i));
         }
      }
   }
   return 0;
}

int get_plane(double *com_buff, double *double_buff, int size_x, int size_y, int plane) {
// Вытаскиваем одну из комплексных плоскостей
// plane - действительная/мнимая (0/1) части комплексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=*(com_buff+2*(j*size_x+i)+plane);
      }
   }
   return 0;
}

int get_plane_abs(double *com_buff, double *double_buff, int size_x, int size_y, int plane) {
// Вытаскиваем одну из комплексных плоскостей
// plane - действительная/мнимая (0/1) части комплексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=fabs(*(com_buff+2*(j*size_x+i)+plane));
//         printf(" %2.4f", *(double_buff+j*size_x+i));         
      }
   }
   return 0;
}

int is_bin(unsigned char *bmp_buff, int size_x, int size_y) {
// Определяем бинарный битмап или нет.
   int i, j;
   int num_col=0;
   unsigned int colors[256];

   for (j=0; j<256; j++) {
      colors[j]=0;
   }

   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         colors[*(bmp_buff+j*size_x+i)]++;
      }
   }

   for (j=0; j<256; j++) {
      if (colors[j]>0) num_col++;
//      printf("%d ", colors[j]);
   }

   if (num_col==2) {return 1;}
   else {return 0;}
}

int make_complex(double *com_buff, double *double_buff_amp, double *double_buff_phs, int size_x, int size_y) {
// Образуем комплексное число через экспоненциальную запись
// (переводим из экспоненциальной записи в координатную)
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(com_buff+2*(j*size_x+i)+PLANE_RE)=(*(double_buff_amp+j*size_x+i))*cos(*(double_buff_phs+j*size_x+i));
         *(com_buff+2*(j*size_x+i)+PLANE_IM)=(*(double_buff_amp+j*size_x+i))*sin(*(double_buff_phs+j*size_x+i));
//         printf(" %2.2f", *(com_buff+2*(j*size_x+i)+PLANE_IM));
      }
   }
   return 0;
}

int norm_2pi_minus_pi(double *double_buff, int size_x, int size_y) {
// Нормализуем фазу по phs=phs*2pi-pi
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=(*(double_buff+j*size_x+i))*2*M_PI-M_PI;
      }
   }
   return 0;
}

int norm_avg(double *double_buff, int size_x, int size_y, double avg1, double avg2) {
// Нормализация по двум средним
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=*(double_buff+j*size_x+i)/avg2;
//         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=*(double_buff+j*size_x+i)*avg1;
//         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int norm_com(double *com_buff, int size_x, int size_y, int plane, double val) {
// Нормализуем комплексный массив по значению val
// plane - действительная/мнимая (0/1) части комплексной матрицы
// plane=2 - нормализация обеих плоскостей сразу
   int i, j;
   double max_re=*(com_buff+2*(0*size_x+0)+PLANE_RE);
   double max_im=*(com_buff+2*(0*size_x+0)+PLANE_IM);
   if (plane==PLANE_BOTH) {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            if (*(com_buff+2*(j*size_x+i)+PLANE_RE) > max_re) {
               max_re=*(com_buff+2*(j*size_x+i)+PLANE_RE);
            }
            if (*(com_buff+2*(j*size_x+i)+PLANE_IM) > max_im) {
               max_im=*(com_buff+2*(j*size_x+i)+PLANE_IM);
            }
         }
      }
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(com_buff+2*(j*size_x+i)+PLANE_RE)=(*(com_buff+2*(j*size_x+i)+PLANE_RE)/max_re)*val;
//            *(com_buff+2*(j*size_x+i)+PLANE_IM)=(*(com_buff+2*(j*size_x+i)+PLANE_IM)/max_im)*val;
         }
      }
   } else {
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            if (*(com_buff+2*(j*size_x+i)+plane) > max_re) {
               max_re=*(com_buff+2*(j*size_x+i)+plane);
            }
         }
      }
      for (j=0; j<size_y; j++) {
         for (i=0; i<size_x; i++) {
            *(com_buff+2*(j*size_x+i)+plane)=(*(com_buff+2*(j*size_x+i)+plane)/max_re)*val;
         }
      }
   }
   return 0;
}

int norm_double(double *double_buff, int size_x, int size_y, double val) {
// Нормализуем вещественный массив по значению val
   int i, j;
   double max=*(double_buff+0*size_x+0);
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         if (*(double_buff+j*size_x+i) > max) {
            max=*(double_buff+j*size_x+i);
         }
      }
   }
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=(*(double_buff+j*size_x+i)/max)*val;
//         printf(" %2.4f", *(double_buff+j*size_x+i));
      }
   }
   return 0;
}

int norm_fft(double *com_buff, int size_x, int size_y) {
// Нормализация матрицы на n=size_x*size_y после fft-1
   int i, j;
   int n=size_x*size_y;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(com_buff+2*(j*size_x+i)+PLANE_RE)=*(com_buff+2*(j*size_x+i)+PLANE_RE)/n;
//         *(com_buff+2*(j*size_x+i)+PLANE_IM)=*(com_buff+2*(j*size_x+i)+PLANE_IM)/n;
         *(com_buff+2*(j*size_x+i)+PLANE_IM)=*(com_buff+2*(j*size_x+i)+PLANE_IM)/n*pow((-1), i*j);
      }
   }
   return 0;
}

int norm_pi_minus_pi(double *double_buff, int size_x, int size_y) {
// Нормализуем фазу по phs=phs*2pi-pi
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(double_buff+j*size_x+i)=(*(double_buff+j*size_x+i))*M_PI-M_PI;
      }
   }
   return 0;
}

int set_plane(double *com_buff, double *double_buff, int size_x, int size_y, int plane) {
// Заполняем одну из комплексных плоскостей
// plane - действительная/мнимая (0/1) части комплексной матрицы
   int i, j;
   for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
         *(com_buff+2*(j*size_x+i)+plane)=*(double_buff+j*size_x+i);
//         printf(" %2.4f", *(double_buff+j*size_x+i));         
      }
   }
   return 0;
}

