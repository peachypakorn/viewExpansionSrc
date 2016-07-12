#include <hls_stream.h>
#include "forwardfftall.h"
#include "nufft_cosfilter.h"
#include "hls_dsp.h"
//#include "hls_cmpy.h"

#include "ap_fixed.h"
#include "hls_fft.h"

// configurable params
const char FFT_INPUT_WIDTH                     = 16;
const char FFT_OUTPUT_WIDTH                    = 32;
const char FFT_CONFIG_WIDTH                    = 8;

#include <complex>
using namespace std;

struct config10 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 16;

    static const unsigned input_width = 16;
    static const unsigned output_width = 32;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 10;
    static const bool has_nfft = true;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};

typedef ap_fixed<FFT_INPUT_WIDTH,1> data_in_t;
typedef ap_fixed<FFT_OUTPUT_WIDTH,FFT_OUTPUT_WIDTH-FFT_INPUT_WIDTH+1> data_out_t;
typedef std::complex<data_in_t> cmpxDataIn;
typedef std::complex<data_out_t> cmpxDataOut;

//#define S 512


template< typename config_f, typename Tin, int FFT_LENGTH>
void dummy_proc_fe(
    bool direction,
	hls::ip_fft::config_t<config_f> * config,
    Tin in[FFT_LENGTH],
    Tin out[FFT_LENGTH],
	int m)
{
#pragma HLS INLINE
#pragma HLS interface ap_fifo port=config
    int i;
    config->setDir(direction);

    if(m == 512)config->setNfft(9);
    else if(m == 256)config->setNfft(8);
    else if(m == 128)config->setNfft(7);
    else if(m == 64)config->setNfft(6);
    else config->setNfft(10);
    for (i=0; i< FFT_LENGTH; i++)
//#pragma HLS INLINE
out[i] = in[i];
}

template< typename config_f, typename Tout, int FFT_LENGTH>
void dummy_proc_be(
		hls::ip_fft::status_t<config_f> * status_in,
    bool* ovflo,
	Tout in[FFT_LENGTH],
	Tout out[FFT_LENGTH])
{

    int i;
    for (i=0; i< FFT_LENGTH; i++) out[i] = in[i];

}

template < typename config_f, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>
void fft_top(
    bool direction,
    Tin in[FFT_LENGTH],
    Tout out[FFT_LENGTH],
    bool* ovflo,
	int length)
{
//#pragma HLS inline

#pragma HLS interface ap_hs port=direction
#pragma HLS interface ap_fifo depth=1 port=ovflo
#pragma HLS interface ap_fifo depth=FFT_LENGTH port=in,out

#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out
#pragma HLS DATAFLOW

	Tin  xn[FFT_LENGTH];
    Tout xk[FFT_LENGTH];
#pragma HLS INTERFACE ap_fifo depth = 1024 port=xn
#pragma HLS INTERFACE ap_fifo depth = 1024 port=xk

#pragma HLS data_pack variable=xn
#pragma HLS data_pack variable=xk
#pragma HLS STREAM variable=xn depth=1024 dim=1
#pragma HLS STREAM variable=xk depth=1024 dim=1
    hls::ip_fft::config_t<config_f> fft_config;

#pragma HLS DATA_PACK variable=fft_config

    hls::ip_fft::status_t<config_f> fft_status;
//#pragma HLS interface ap_fifo port=fft_config.data

    int m;
    if(length==1024)m=1024;
    else if(length==512)	m=512;
    else if(length==256)m = 256; //dummy_proc_fe<config_f, Tin, 256,FFT_LENGTH>(direction, &fft_config, in, xn);
    else if(length==128)m = 128; //dummy_proc_fe<config_f, Tin, 128,FFT_LENGTH>(direction, &fft_config, in, xn);
    else if(length==64)m  = 64;				 //dummy_proc_fe<config_f, Tin, 64,FFT_LENGTH>(direction, &fft_config, in, xn);
    dummy_proc_fe<config_f, Tin,FFT_LENGTH>(direction, &fft_config, in, xn,m);
    // FFT IP
    hls::fft<config_f>(xn, xk, &fft_status, &fft_config);
    dummy_proc_be<config_f, Tout, FFT_LENGTH>(&fft_status, ovflo, xk, out);

}






template < typename config_f, int DEPTHIN, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>
void fft_top_stream1( hls::stream< t_nufft_output_complex > &in,
					hls::stream< t_nufft_output_complex > &outH,
					hls::stream< t_nufft_output_complex > &outL,
					 hls::stream< t_nufft_output_complex > &out512,
					 hls::stream< t_nufft_output_complex > &out256,
					 hls::stream< t_nufft_output_complex > &out128,
					 hls::stream< t_nufft_output_complex > &out64
					 ) {
//#pragma HLS INLINE
#pragma HLS data_pack variable=in
#pragma HLS interface ap_fifo depth=1024 port=outH
#pragma HLS interface ap_fifo depth=1024 port=outL
#pragma HLS interface ap_fifo depth=512 port=out512
#pragma HLS interface ap_fifo depth=256 port=out256
#pragma HLS interface ap_fifo depth=128 port=out128
#pragma HLS interface ap_fifo depth=64 port=out64

//#pragma HLS data_pack variable=out512
//#pragma HLS data_pack variable=out256
//#pragma HLS data_pack variable=out128
//#pragma HLS data_pack variable=out64

#pragma HLS DATAFLOW

	Tin  inM[FFT_LENGTH];
	Tout outM[FFT_LENGTH];
#pragma HLS interface ap_fifo depth=1024 port=inM,outM
#pragma HLS data_pack variable=inM
#pragma HLS data_pack variable=outM
//#pragma HLS STREAM variable=DEPTHIN
	t_nufft_output_complex tempcheck[1024+1024+960];//!!!!!!!!!!!Problem

	int length[6]={1024,1024,512,256,128,64};
	for (int level = 0; level < 6; level++) {
	for(int i=0;i<length[level];i++) {
#pragma HLS PIPELINE
		t_nufft_output_complex valIn;
		if(i<length[level]){
			valIn = in.read();
			inM[i] = cmpxDataIn(valIn.real(), valIn.imag());
			}
		}

	bool direction;
	bool ovflo;
	//const int LEN = FFT_LENGTH>>level;

	fft_top<config_f,FFT_LENGTH, Tin, Tout, unrelated_manner>(true, inM, outM, &ovflo,length[level]);

	int len[6] = {0,1024,1024+1024, 2048+512,2048+512+256,2048+512+256+128};
	for(int i=0;i<FFT_LENGTH;i++) {
#pragma HLS PIPELINE
		Tout valOut = outM[i];

		t_nufft_output_complex val(valOut.real(), valOut.imag());
		 tempcheck[i+len[level]] = val;
		 //float a = val.real().to_float();
		if(level==0&&i<length[level])outH.write(val);
		else if(level==1&&i<length[level])outL.write(val);
		else if(level==2&&i<length[level])out512.write(val);
		else if(level==3&&i<length[level])out256.write(val);
		else if(level==4&&i<length[level])out128.write(val);
		else if(level==5&&i<length[level])out64.write(val);
	}
	}
#ifndef __SYNTHESIS__
	{

		FILE * fo = fopen("fftOut.txt", "wb");
		for(int i=0;i<3008;i++) {
			fprintf(fo, "%.8f %.8f\n",  tempcheck[i].real().to_double(),  tempcheck[i].imag().to_float());
		}
		fclose(fo);

	}
#endif
}


void s2_fft_all_stream(
		hls::stream< t_nufft_output_complex > &inH,
		hls::stream< t_nufft_output_complex > &inL0,
		hls::stream< t_nufft_output_complex > &inLA,
		hls::stream< t_nufft_output_complex > &outH,
		hls::stream< t_nufft_output_complex > &outL0,
		hls::stream< t_nufft_output_complex > &outLA
		){
#pragma HLS INTERFACE axis port=inH
#pragma HLS INTERFACE axis port=inL0
#pragma HLS INTERFACE axis port=inLA
#pragma HLS INTERFACE axis port=outH
#pragma HLS INTERFACE axis port=outL0
#pragma HLS INTERFACE axis port=outLA
#pragma HLS data_pack variable=inH
#pragma HLS data_pack variable=inL0
#pragma HLS data_pack variable=inLA
#pragma HLS data_pack variable=outH
#pragma HLS data_pack variable=outL0
#pragma HLS data_pack variable=outLA

#pragma HLS DATAFLOW

	// Assuming that the buffer is actually on the other side
	hls::stream< t_nufft_output_complex > FFTHOut;
	hls::stream< t_nufft_output_complex > FFTL0Out;
		hls::stream< t_nufft_output_complex > FFT512Out;
		hls::stream< t_nufft_output_complex > FFT256Out;
		hls::stream< t_nufft_output_complex > FFT128Out;
		hls::stream< t_nufft_output_complex > FFT64Out;
	#pragma HLS DATA_PACK variable=FFTHOut
	#pragma HLS DATA_PACK variable=FFTL0Out
	#pragma HLS DATA_PACK variable=FFT512Out
	#pragma HLS DATA_PACK variable=FFT256Out
	#pragma HLS DATA_PACK variable=FFT128Out
	#pragma HLS DATA_PACK variable=FFT64Out
	#pragma HLS STREAM variable=FFTHOut depth=1024
	#pragma HLS STREAM variable=FFTL0Out depth=1024
	#pragma HLS STREAM variable=FFT512Out depth=512
	#pragma HLS STREAM variable=FFT256Out depth=256
	#pragma HLS STREAM variable=FFT128Out depth=128
	#pragma HLS STREAM variable=FFT64Out depth=64

		//hls::stream< t_nufft_output_complex > nFFTHOut;
		//hls::stream< t_nufft_output_complex > nFFTL0Out;
		hls::stream< t_nufft_output_complex > nFFT256Out;
		hls::stream< t_nufft_output_complex > nFFT128Out;
		hls::stream< t_nufft_output_complex > nFFT64Out;
		hls::stream< t_nufft_output_complex > nFFT32Out;

		//#pragma HLS DATA_PACK variable=nFFTHOut
		//#pragma HLS DATA_PACK variable=nFFTL0Out
		#pragma HLS DATA_PACK variable=nFFT256Out
		#pragma HLS DATA_PACK variable=nFFT128Out
		#pragma HLS DATA_PACK variable=nFFT64Out
		#pragma HLS DATA_PACK variable=nFFT32Out

	//#pragma HLS STREAM variable=nFFTHOut  depth=512
	//#pragma HLS STREAM variable=nFFTL0Out  depth=512
	#pragma HLS STREAM variable=nFFT256Out  depth=256
	#pragma HLS STREAM variable=nFFT128Out  depth=128
	#pragma HLS STREAM variable=nFFT64Out  depth=64
	#pragma HLS STREAM variable=nFFT32Out  depth=32
		hls::stream< t_nufft_output_complex > sigIn;
#pragma HLS DATA_PACK variable=sigIn
#pragma HLS STREAM variable=sigIn depth=3008 dim=1

for (int i = 0; i < 3008; i++) {
	t_nufft_output_complex tmp;
	if(i<1024){tmp = inH.read();
		tmp.real() = tmp.real()*4;
		tmp.imag() = tmp.imag()*4;
	}
	else if(i>=1024&&i<2048){
		tmp = inL0.read();
		tmp.real() = tmp.real()*4;
		tmp.imag() = tmp.imag()*4;
	}
	else if(i>=2048&&i<3008) {
		tmp = inLA.read();
		tmp.real() = tmp.real()/2;
		tmp.imag() = tmp.imag()/2;
	}
		sigIn.write(tmp);
}

		fft_top_stream1< config10, 1536, 1024, cmpxDataIn, cmpxDataOut, 0> ( sigIn,FFTHOut,FFTL0Out, FFT512Out,FFT256Out,FFT128Out,FFT64Out);

				nufft_cosfilter1024<3, t_nufft_output_complex>(FFTHOut,outH);
				nufft_cosfilter1024<3, t_nufft_output_complex>(FFTL0Out, outL0);
				nufft_cosfilter_single<t_nufft_output_complex>(FFT512Out, nFFT256Out, 256, 2);
				nufft_cosfilter_single<t_nufft_output_complex>(FFT256Out, nFFT128Out, 128, 2);
				nufft_cosfilter_single<t_nufft_output_complex>(FFT128Out, nFFT64Out, 64, 2);
				nufft_cosfilter_single<t_nufft_output_complex>(FFT64Out, nFFT32Out, 32, 2);

				LoopDataOut:
				for(int i=0;i<(256+128+64+32);i++) {
			#pragma HLS pipeline
					t_nufft_output_complex  v0;
					if (i<256) {
						v0 = nFFT256Out.read();

					}
					else if (i>=256 && i < 384) {
						v0 = nFFT128Out.read();

					}
					else if (i>=384 && i< 448) {
						v0 = nFFT64Out.read();

					}
					else if (i>=448&& i< 480) {
						v0 = nFFT32Out.read();

					}
					outLA.write(v0);

				}

}



