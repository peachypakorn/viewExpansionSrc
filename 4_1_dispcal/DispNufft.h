#ifndef __DISPNUFFT__
#define __DISPNUFFT__

#include "common_type.h"
//MedianFiler
#include "../3_MedianFilter/median_filter.h"
#include "ap_fixed.h"
#include "hls_fft.h"
#include <hls_stream.h>
//#include <hls_video.h>
#include <ap_axi_sdata.h>

void nufft_Disp (
		hls::stream<data_t>					&src_axi,
		hls::stream<t_disp_scalar>          &dispFilters);

#endif
