//---------------------------------------------------------------------------

#ifndef fncH
#define fncH
//---------------------------------------------------------------------------
#endif

int bitmap_to_double(double *double_buff, unsigned char *bmp_buff, int size_x, int size_y);
int com_bitmap_to_plane(double *com_buff, unsigned char *bmp_buff, int size_x, int size_y, int plane);
int com_copy_image(double *com_src, double *com_dst, int size_x, int size_y);
int com_echo(double *com_buff, int size_x, int size_y);
int com_fill_plane(double *com_buff, int size_x, int size_y, int plane, int val, int rnd);
int com_get_amp(double *com_buff, double *double_buff, int size_x, int size_y);
int com_get_phase(double *com_buff, double *double_buff, int size_x, int size_y);
int com_plane_to_bitmap(double *com_buff, unsigned char *bmp_buff, int size_x, int size_y, int plane);
double com_rmse(double *double_buff1, double *double_buff2, int size_x, int size_y);
int com_shift(double *com_buff, int size_x, int size_y);
double double_avg(double *double_buff, int size_x, int size_y, int absol);
int double_binarize(double *double_buff, int size_x, int size_y);
int double_echo(double *double_buff, int size_x, int size_y);
int double_fill_plane(double *double_buff, int size_x, int size_y, int val, int rnd);
int double_rescale(double *double_buff, int size_x, int size_y, int c_min, int c_max);
int double_shift(double *double_buff, int size_x, int size_y);
int double_sinc(double *double_buff, int size_x, int size_y, double corr);
int double_to_bitmap(double *double_buff, unsigned char *bmp_buff, int size_x, int size_y);
int fill_double(double *double_buff, int size_x, int size_y, int val, int rnd);
int get_plane(double *com_buff, double *double_buff, int size_x, int size_y, int plane);
int get_plane_abs(double *com_buff, double *double_buff, int size_x, int size_y, int plane);
int is_bin(unsigned char *bmp_buff, int size_x, int size_y);
int make_complex(double *com_buff, double *double_buff_amp, double *double_buff_phs, int size_x, int size_y);
int norm_2pi_minus_pi(double *double_buff, int size_x, int size_y);
int norm_avg(double *com_buff, int size_x, int size_y, double avg1, double avg2);
int norm_com(double *com_buff, int size_x, int size_y, int plane, double val);
int norm_double(double *double_buff, int size_x, int size_y, double val);
int norm_fft(double *com_buff, int size_x, int size_y);
int norm_pi_minus_pi(double *double_buff, int size_x, int size_y);
int set_plane(double *com_buff, double *double_buff, int size_x, int size_y, int plane);
