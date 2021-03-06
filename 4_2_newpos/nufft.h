#ifndef __NUFFT___
#define __NUFFT___

//#include "../DispCheck/cmpy_complex.h"
//#include "../DispCheck/sincos.h"
#include "common_type.h"

//typedef ap_fixed<17, 4> t_nufft_output_scalar;
//typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;

void nufft_top_pyr(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOutH,
			   hls::stream<t_nufft_output_complex> &sigStreamOutL0,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLA,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLP) ;

void nufft_interface(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOutH,
			   hls::stream<t_nufft_output_complex> &sigStreamOutL0,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLA,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLP);

#endif
