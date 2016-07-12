#include "hls_dsp.h"
#include <fstream>
#include <complex>
#include "hls_fft.h"

#include "pyr.h"

#if 1
const int Kset = 8;
const int limits[] = { 512,512, 512,256,128,64,32,16};
const int climits[] = { 0,0, 512, 1024,1280,1408,1472,1504};
const int llimits[] = { 9,9,9,8,7,6,5,4};
const int lshifts[] = { 0,1, 2, 3,3,4,5,9};
const int pshifts[] = { 0,8, 7, 5,4,2,0, 0};
const int fsizes[]={ 512,512,512,256,128,64,32,16};
#endif

template <typename T>
void pyrbuild_top(cmpxDataOut  fftTmp[512],
				  T   fftPyrFilOut[ 1520],
			      int width,
                  const int nL){
	#pragma HLS INTERFACE ap_fifo depth=512 port=fftTmp
	#pragma HLS INTERFACE ap_fifo depth=1520 port=fftPyrFilOut
//#pragma HLS STREAM variable=fftTmp depth=2048 dim=1
//#pragma HLS STREAM variable=fftPyrFilOut depth=1024 dim=1
#pragma HLS INLINE
//#pragma HLS INTERFACE axis port=fftIn
//#pragma HLS INTERFACE axis port=pyrFilOut
#if 0

#else

	int l = 0;
	int i = 0;
	for(int coefIdx = 0;coefIdx < 1520; coefIdx++)
	{
#pragma HLS PIPELINE
		int nlimit = limits[l+1];
		int half = nlimit >> 1;
		int idx;
		T val ;
		if (i<half) idx = i;
		//i also have value
		else  idx = 512- nlimit + i;
		val = fftTmp[idx] ;
		T outVal;

//		std::cout << i << " "<< idx << "  " << coefIdx << std::endl;

		outVal.real() = val.real() * consFilters[coefIdx];
		outVal.imag() = val.imag() * consFilters[coefIdx];
		fftPyrFilOut[coefIdx] = outVal;
		//printf("%d %d %d %f %f \n",coefIdx,i,idx,outVal.real().to_float(),consFilters[coefIdx].to_float());
		i++;
		l += (nlimit == i);
		i = (nlimit == i)?0:i;
	}

#endif
}
#define IMG_WIDTH    512


template< typename config_f, typename Tin, int FFT_LENGTH>
void dummy_proc_fe(
    bool direction,
	hls::ip_fft::config_t<config_f> * config,
    Tin in[FFT_LENGTH],
    Tin out[FFT_LENGTH],
	int m)
{
	//#pragma HLS INTERFACE ap_fifo depth=1024 port=in
	//#pragma HLS INTERFACE ap_fifo depth=1024 port=out
//#pragma HLS INLINE
#pragma HLS interface ap_fifo port=config
    int i;
    config->setDir(direction);
    if(m == 256)config->setNfft(8);
    else if(m == 128)config->setNfft(7);
    else if(m == 64)config->setNfft(6);
    else if(m == 32)config->setNfft(5);
    else if(m == 16)config->setNfft(4);
    else config->setNfft(9);
   // config->setNfft(log_size);
//    config->setSch(0x2AB);
    for (i=0; i< FFT_LENGTH; i++)
#pragma HLS PIPELINE
out[i] = in[i];
}


template< typename T>
void dummy_proc2( T imgIn[IMG_WIDTH], T out[IMG_WIDTH]){
#pragma HLS inline
	for(int i=0;i<512;i++)
		out[i] = imgIn[i];
}

template< typename T>
void dummy_proc3( T imgIn[IMG_WIDTH],hls::stream<T> &out,int length){
#pragma HLS inline
	for(int i=0;i<length;i++)
		out.write( imgIn[i]);
}

template< typename config_f, typename Tout, int FFT_LENGTH>
void dummy_proc_be(
		hls::ip_fft::status_t<config_f> * status_in,
    bool* ovflo,
	Tout in[FFT_LENGTH],
	Tout out[FFT_LENGTH])
{

    int i;
    for (i=0; i< FFT_LENGTH; i++)out[i] = in[i];
//    *ovflo = status_in->getOvflo() & 0x1;
}


template < typename config_f, int FFT_LENGTH, typename Tin, typename Tout>
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
#pragma HLS interface ap_fifo depth=512 port=in,out

#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out
#pragma HLS DATAFLOW
	Tin  xn[FFT_LENGTH];
    Tout xk[FFT_LENGTH];
#pragma HLS INTERFACE ap_fifo depth = 512 port=xn
#pragma HLS INTERFACE ap_fifo depth = 512 port=xk
#pragma HLS data_pack variable=xn
#pragma HLS data_pack variable=xk
#pragma HLS STREAM variable=xn depth=512 dim=1
#pragma HLS STREAM variable=xk depth=512 dim=1
    hls::ip_fft::config_t<config_f> fft_config;
//#pragma HLS FUNCTION_INSTANTIATE variable=fft_config

#pragma HLS DATA_PACK variable=fft_config
//#pragma HLS DATA_PACK variable=fft_config
    hls::ip_fft::status_t<config_f> fft_status;
//#pragma HLS interface ap_fifo port=fft_config
    int m;
    if(length==512)	m=512;
    else if(length==256)m = 256; //dummy_proc_fe<config_f, Tin, 256,FFT_LENGTH>(direction, &fft_config, in, xn);
    else if(length==128)m = 128; //dummy_proc_fe<config_f, Tin, 128,FFT_LENGTH>(direction, &fft_config, in, xn);
    else if(length==64)m = 64;
    else if(length==32)m = 32;
    else m  = 16;				 //dummy_proc_fe<config_f, Tin, 64,FFT_LENGTH>(direction, &fft_config, in, xn);
    dummy_proc_fe<config_f, Tin, FFT_LENGTH>(direction, &fft_config, in, xn, m);
    // FFT IP
    hls::fft<config_f>(xn, xk, &fft_status, &fft_config);
    dummy_proc_be<config_f, Tout, 512>(&fft_status, ovflo, xk, out);

}

void pyrconstuct_top(
			  cmpxDataIn imgIn[IMG_WIDTH],
		      //hls::stream<t_image> &imgIn,
			  hls::stream<t_pyr_complex> &pyrFilOut,
		      const int nL
		      ) {
#pragma HLS DATAFLOW

#pragma HLS interface ap_fifo depth=512 port=imgIn
#pragma HLS interface ap_fifo depth=1520 port=pyrFilOut
#pragma HLS data_pack variable=pyrFilOut
#pragma HLS data_pack variable=imgIn
	hls::stream<cmpxDataIn> imgInP;
	cmpxDataIn imgInTmp[FFT_LENGTH];
	cmpxDataOut imgOutTmpFFTStream[FFT_LENGTH];
	cmpxDataOut imgOutTmpBlockRam[FFT_LENGTH];
	cmpxDataOut  fftPyrOut[1520];
	cmpxDataOut LP[16];

//#pragma HLS STREAM variable=imgOutTmpFFTStream depth=1024 dim=1
//#pragma HLS STREAM variable=fftPyrOut depth=2048 dim=1
#pragma HLS interface ap_fifo depth=512 port=imgOutTmpFFTStream
#pragma HLS interface ap_fifo depth=512 port=imgOutTmpBlockRam
#pragma HLS interface ap_fifo depth=512 port=imgInTmp
#pragma HLS interface ap_fifo depth=1520 port=fftPyrOut
	#pragma HLS interface ap_fifo depth=512 port=imgInP
	#pragma HLS data_pack variable=imgInP
	#pragma HLS STREAM variable=imgInP depth=512 dim=1
//	for (int i = 0; i < 512; i++) {
//		imgInP[i] = imgIn[i];
//	}
	dummy_proc3(imgIn,imgInP,512);
	LPH: for(int l=0;l<8;l++){
//#pragma HLS DATAFLOW
		const int	cidx = climits[l];
		const int	nlimit = limits[l];
		const int	fsize = fsizes[l];
		const int	lshift = lshifts[l];
		bool direction =false;
		bool ovflo;
		if( l ==0){
					direction = true;
					//for(int i=0;i<512;i++) {
					//#pragma HLS pipeline
								//cmpxDataIn val (0,0);
									//val.real() = imgIn[i].real();
									//val.imag() =imgIn[i].real()
								//imgInTmp[i] = imgIn[i];
							//}
				  }

			for(int i=0;i<512;i++) {
			#pragma HLS pipeline
						cmpxDataIn val (0,0);
						if(l==0){
							//val = imgInP[i];
							val = imgInP.read();
							//val.real() = imgIn[i].real();
							//val.imag() = imgIn[i].imag();
						}
						else if (i<nlimit) {
							val.real() = (fftPyrOut[cidx + i ].real() >> lshift);
							val.imag() = (fftPyrOut[cidx + i ].imag() >> lshift);
						}
						if(l==7&&i<nlimit) LP[i] = val;
						imgInTmp[i] = val;
					}
#ifndef __SYNTHESIS__
	if (l==1)
	{
		FILE * fo = fopen("ifft_in.txt", "wb");
		for(int i=0;i<512;i++) {
			fprintf(fo, "%.8f %.8f\n", imgInTmp[i].real().to_float(), imgInTmp[i].imag().to_float());
		}
		fclose(fo);


	}
#endif

		fft_top<config2,IMG_WIDTH,cmpxDataIn,cmpxDataOut>(direction,imgInTmp,imgOutTmpFFTStream,&ovflo,nlimit);
		dummy_proc2(imgOutTmpFFTStream, imgOutTmpBlockRam);
		if(l==0){

			pyrbuild_top(imgOutTmpBlockRam, fftPyrOut, 512, 512);

#ifndef __SYNTHESIS__
				{
					FILE * fo = fopen("fftOut.txt", "wb");
					for(int i=0;i<512;i++) {
						fprintf(fo, "%.8f      %.8f\n", imgOutTmpFFTStream[i].real().to_float(), imgOutTmpFFTStream[i].imag().to_float());
					}
					fclose(fo);
				}
				{
					FILE * fo = fopen("fftOutFilter.txt", "wb");
					for(int i=0;i<1520;i++) {
						fprintf(fo, "%.8f       %.8f\n", fftPyrOut[i].real().to_float(), fftPyrOut[i].imag().to_float());
					}
				}
			/*
				for(int i=0;i<512;i++)
					std::cout << "FFT Out " << i << " : " << imgOutTmp[i] << std::endl;

			*/
#endif
		}//end of input
		else{
			cmpxDataOut first = imgOutTmpFFTStream[0];
		for(int i=1;i<=512;i++) {
#pragma HLS pipeline
			if(i<=nlimit){
			t_pyr_complex val;
			cmpxDataOut tmp;
			if(i==nlimit)tmp = first;
			else tmp = imgOutTmpFFTStream[i];
			val.real() = tmp.real() >> pshifts[l];
			val.imag() = tmp.imag() >> pshifts[l];
			if(l==7)val = cmpxDataOut(LP[i-1].real(),LP[i-1].imag());
			pyrFilOut.write(val);
				}
			}
		}
		//end of ifft
	}
}

void pyrcon_top(
				hls::stream<cmpxDataIn> &imgIn,
		      hls::stream<t_pyr_complex> &pyrFilOut,
		      const int nL
		      ){
#pragma HLS STREAM variable=imgIn depth=512 dim=1
#pragma HLS DATA_PACK variable=pyrFilOut
#pragma HLS DATA_PACK variable=imgIn
#pragma HLS DATAFLOW
#pragma HLS INTERFACE axis depth=1520 port=pyrFilOut
#pragma HLS INTERFACE axis depth=512 port=imgIn

	cmpxDataIn In0[IMG_WIDTH];
#pragma HLS STREAM variable=In0 depth=512 dim=1
#pragma HLS INTERFACE ap_fifo depth=512 port=In0
	hls::stream<t_pyr_complex> Out;
#pragma HLS STREAM variable=Out depth=1520 dim=1
#pragma HLS INTERFACE ap_fifo depth=512 port=Out
	int i ;
	for (i = 0; i < 512; i++) {
#pragma HLS PIPELINE
		In0[i] = imgIn.read();
	}
	pyrconstuct_top(In0,Out,nL);

	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		t_pyr_complex val = Out.read(); ;
		pyrFilOut.write(val);
	}


}


