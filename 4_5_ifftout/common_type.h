#ifndef __COMMON_TYPE
#define __COMMON_TYPE

#include "ap_fixed.h"
#include "hls_fft.h"
#include <hls_stream.h>
//#include <hls_video.h>
#include <ap_axi_sdata.h>


#include <complex>

typedef ap_fixed<17, 4> t_nufft_output_scalar;
typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;

typedef ap_axiu<112, 1, 1, 1> PACKETIN;
typedef ap_axiu<24, 1, 1, 1> PIXEL;


typedef  ap_uint<112>  packed;

typedef ap_uint<24> PIXEL_RAW;

typedef hls::stream<PACKETIN> AXI_STREAMIN;

typedef hls::stream<PIXEL_RAW > AXI_STREAMRGB;

//typedef hls::Scalar<3, unsigned char> RGB_PIXEL;
//Configuration parameters for this instance.
//typedef hls::CmpyThreeMult ARCH;
//ROUND_MODE Values : AP_RND, AP_RND_ZERO, AP_RND_MIN_INF, AP_RND_INF, AP_RND_CONV, AP_TRN, AP_TRN_ZERO
//OVERFLOW_MODE Values : AP_SAT, AP_SAT_ZERO, AP_SAT_SYM, AP_WRAP, AP_WRAP_SM

const int INPUT_WIDTH                = 20;
const int INPUT_INTEGER_BITS         = 3;
const int INPUT_ROUND_MODE           = 5;//AP_TRN, etc;
const int INPUT_OVERFLOW_MODE        = 3;//AP_WRAP, etc;
const int INPUT_SATURATION_BITS      = 0;

const int OUTPUT_WIDTH               = 20;
const int OUTPUT_INTEGER_BITS        = 2;
const int OUTPUT_ROUND_MODE          = 0;//AP_TRN, etc;
const int OUTPUT_OVERFLOW_MODE       = 3;//AP_WRAP, etc;
const int OUTPUT_SATURATION_BITS     = 0;

//Interface types
typedef ap_fixed<
  INPUT_WIDTH,
  INPUT_INTEGER_BITS,
  (ap_q_mode)INPUT_ROUND_MODE,
  (ap_o_mode)INPUT_OVERFLOW_MODE,
  INPUT_SATURATION_BITS> t_input_scalar;


const int DISP_WIDTH                = 24;
const int DISP_INTEGER_BITS         = 12;
const int DISP_ROUND_MODE           = 5;//AP_TRN, etc;
const int DISP_OVERFLOW_MODE        = 3;//AP_WRAP, etc;
const int DISP_SATURATION_BITS      = 0;

//Interface types
typedef ap_fixed<
		DISP_WIDTH,
		DISP_INTEGER_BITS,
  (ap_q_mode)DISP_ROUND_MODE,
  (ap_o_mode)DISP_OVERFLOW_MODE,
  DISP_SATURATION_BITS> t_disp_scalar;


typedef ap_fixed<
  OUTPUT_WIDTH,
  OUTPUT_INTEGER_BITS,
  (ap_q_mode)OUTPUT_ROUND_MODE,
  (ap_o_mode)OUTPUT_OVERFLOW_MODE,
  OUTPUT_SATURATION_BITS> t_output_scalar;

typedef std::complex< t_input_scalar > t_input_complex;
typedef std::complex< t_output_scalar > t_output_complex;

typedef ap_fixed<15, 3>   t_coeff_scalar;
typedef std::complex < t_coeff_scalar> t_coeff;



inline t_nufft_output_complex fromAP( ap_uint<30> v) {
	t_coeff_scalar r, i;

	r.range() = v.range(29,15);
	i.range() = v.range(14,0);

	t_nufft_output_complex vOut( r, i);

	return vOut;
}



typedef ap_fixed<24,10>            t_recon_scalar;
typedef std::complex< t_recon_scalar>   t_recon_complex;






#ifndef __SYNTHESIS__
template <typename T>
void DataSniffer(hls::stream<T> &t, std::vector<T> &vecOut) {
	vecOut.clear();
	while(!t.empty())
		vecOut.push_back(t.read());
	for(int i=0;i<(int)vecOut.size();i++) t.write(vecOut[i]);
}


template <typename T>
void DataPlot(const char *name, std::vector<T> &vec) {
	FILE * fo = fopen(name, "w");
	for(int i=0;i<(int)vec.size();i++)
		MyOut(fo, vec[i]);
	fclose(fo);

}

template <typename T>
void WritePlot(const char *name, hls::stream<T> &t) {
	std::vector<T> vec;
	DataSniffer(t, vec);
	DataPlot( name, vec);

}
#endif

#endif
