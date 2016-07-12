#ifndef __NUFFT___
#define __NUFFT___

//#include "../DispCheck/cmpy_complex.h"
//#include "../DispCheck/sincos.h"
#include "common_type.h"

//typedef ap_fixed<17, 4> t_nufft_output_scalar;
//typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;

void nufft_top_pyr(hls::stream<t_input_complex>  sig[3],
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> sigStreamOutH[3],
			   hls::stream<t_nufft_output_complex> sigStreamOutL0[3],
			   hls::stream<t_nufft_output_complex> sigStreamOutLA[3],
			   hls::stream<t_nufft_output_complex> sigStreamOutLP[3]) ;

#endif
