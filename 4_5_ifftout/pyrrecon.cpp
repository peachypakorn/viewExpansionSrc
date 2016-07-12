#include "common_type.h"
#include <hls_fft.h>
#include "pyrrecon.h"

//#include <hls_video.h>

//#define DEBUGOUT

// NUFFT Format
/*
struct config1 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 16;

    static const unsigned input_width = 16;
    static const unsigned output_width = 16;

    static const unsigned scaling_opt = hls::ip_fft::scaled;

    static const unsigned max_nfft = 10;
    static const bool has_nfft = true;
};
*/

//#define C   1
#define MEMSIZE      512


const int Kset = 7;
const int limits[] = { 512, 512,256,128,64,32,16};

void pyr_recon(hls::stream< cmpxDataOut>  &pyrFFT,
		       hls::stream< t_recon_complex>  &imDFTOut) {
#pragma HLS DATAFLOW

#pragma HLS data_pack variable=pyrFFT
#pragma HLS data_pack variable=imDFTOut

	#if 0
#pragma HLS dataflow

	hls::stream< t_recon_complex > strRead0[C];
	hls::stream< t_recon_complex > strRead1[C];
	hls::stream< t_recon_complex > strRead2[C];
	hls::stream< t_recon_complex > strRead3[C];
	hls::stream< t_recon_complex > strRead4[C];
	hls::stream< t_recon_complex > strRead5[C];
	hls::stream< t_recon_complex > strRead6[C];
#pragma HLS STREAM variable=strRead0 depth=512 dim=1
#pragma HLS STREAM variable=strRead1 depth=512 dim=1
#pragma HLS STREAM variable=strRead2 depth=256 dim=1
#pragma HLS STREAM variable=strRead3 depth=128 dim=1
#pragma HLS STREAM variable=strRead4 depth=64 dim=1
#pragma HLS STREAM variable=strRead5 depth=32 dim=1
#pragma HLS STREAM variable=strRead6 depth=16 dim=1

#pragma HLS data_pack variable=strRead0
#pragma HLS data_pack variable=strRead1
#pragma HLS data_pack variable=strRead2
#pragma HLS data_pack variable=strRead3
#pragma HLS data_pack variable=strRead4
#pragma HLS data_pack variable=strRead5
#pragma HLS data_pack variable=strRead6
//#pragma HLS STREAM variable=strRead[0][1] depth=512 dim=1
//#pragma HLS STREAM variable=strRead[0][2] depth=256 dim=1

	int l = 0;
	int i = 0;
	for(int coefIdx = 0;coefIdx < 1520; coefIdx++)
	{
#pragma HLS PIPELINE
		int nlimit = limits[l];
		int half = nlimit >> 1;
		int idx;

		data_out_t   coef = reconsFilters[coefIdx];
		for(int c=0;c<C;c++) {
#pragma HLS unroll
			cmpxDataOut  val = pyrFFT[c].read();
			t_recon_complex  valMul;
			valMul.real()= val.real() * coef;
			valMul.imag()= val.imag() * coef;

			if (l==0) strRead0[c].write(valMul);
			if (l==1) strRead1[c].write(valMul);
			if (l==2) strRead2[c].write(valMul);
			if (l==3) strRead3[c].write(valMul);
			if (l==4) strRead4[c].write(valMul);
			if (l==5) strRead5[c].write(valMul);
			if (l==6) strRead6[c].write(valMul);
			//strRead[c][l].write(valMul);
		}
		i++;
		l += (nlimit == i);
		i = (nlimit == i)?0:i;
	}

	for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
		for(int c=0;c<C;c++) {
#pragma HLS unroll
			t_recon_complex val0 = strRead0[c].read();
			t_recon_complex val1 = strRead1[c].read();
			t_recon_complex val2 = t_recon_complex(0,0);
			t_recon_complex val3 = t_recon_complex(0,0);
			t_recon_complex val4 = t_recon_complex(0,0);
			t_recon_complex val5 = t_recon_complex(0,0);
			t_recon_complex val6 = t_recon_complex(0,0);

			if (i<128 || i>=(512-128)) val2 = strRead2[c].read();
			if (i<64  || i>=(512- 64)) val3 = strRead3[c].read();
			if (i<32  || i>=(512- 32)) val4 = strRead4[c].read();
			if (i<16  || i>=(512- 16)) val5 = strRead5[c].read();
			if (i<8  || i>=(512- 8))   val6 = strRead6[c].read();

			t_recon_complex  outVal = val0 + val1 + val2 + val3 + val4 + val5 + val6;
			imDFTOut[c].write(outVal);
		}
	}


#else
//	t_output_scalar mem[C][512];
	t_recon_complex mem[MEMSIZE];
#pragma HLS ARRAY_PARTITION variable=mem complete factor=3 dim=1
#pragma HLS data_pack variable=mem

#pragma HLS RESOURCE variable=mem core=RAM_2P_BRAM

	int l = 0;
	int i = 0;
	for(int coefIdx = 0;coefIdx < 1520; coefIdx++)
	{
#pragma HLS INLINE
//#pragma HLS PIPELINE
		int nlimit = limits[l];
		int half = nlimit >> 1;
		int idx;
		data_out_t   coef = reconsFilters[coefIdx];
		if (i<half) idx = i;
		else  idx = MEMSIZE- nlimit + i;

		// Might be better to change the coefficient
		if (coefIdx >= 1504) {
			coef = coef << 6;
		}

		//for(int c=0;c<reconC;c++) {
//#pragma HLS unroll

			cmpxDataOut  val = pyrFFT.read();
			t_recon_complex  valMul;
			valMul.real()= val.real() * coef;
			valMul.imag()= val.imag() * coef;

			/*t_recon_complex  valMem = mem[c][idx];
			t_recon_complex  valSum;
			if (coefIdx < MEMSIZE)
				valSum = valMul + valMem;
			else
				valSum = valMul;
			mem[c][idx] = valSum;*/
			//printf("%d %d %f %f\n", coefIdx, i, val.real().to_float(), coef.to_float());

//			mem[c][idx] = valMul;

			if (coefIdx < MEMSIZE) {
				mem[idx] = valMul;
			} else
				mem[idx] += valMul;

		//}
		i++;
		l += (nlimit == i);
		i = (nlimit == i)?0:i;
	}
	for(int i=0;i<MEMSIZE;i++) {
#pragma HLS PIPELINE
	//	for(int c=0;c<reconC;c++) {
//#pragma HLS unroll
			imDFTOut.write(mem[i]);
		//}
	}
#endif

}


//typedef ap_fixed<16,1> data_in_t;
//typedef ap_fixed<24, 24 - 9 + 1 > data_out_t;

//typedef std::complex<data_in_t> cmpxDataIn;
//typedef std::complex<data_out_t> cmpxDataOut;

void cmifft(  t_recon_complex  DFT[MEMSIZE],
		     ifft_out_t  sigOut[MEMSIZE]) {

#pragma HLS data_pack variable=DFT
#pragma HLS data_pack variable=sigOut
	cmpxIFFTIn   DFTOutMem[MEMSIZE];
	cmpxIFFTOut   IMOutMem[MEMSIZE];

#pragma HLS data_pack variable=DFTOutMem
#pragma HLS data_pack variable=IMOutMem

#pragma HLS STREAM variable=DFTOutMem depth=512 dim=1
#pragma HLS STREAM variable=IMOutMem  depth=512 dim=1

	config2_recon_t  fft_config2;
	status2_recon_t  fft_status2;

#pragma HLS DATAFLOW
	for(int i=0;i<MEMSIZE;i++) {
#pragma HLS PIPELINE
		t_recon_complex t = DFT[i];
		t.real() = t.real()/2;
		t.imag() = t.imag()/2;
		cmpxIFFTIn    val;
		val.real().range(23,0) = t.real().range(23, 2);
		val.imag().range(23,0) = t.imag().range(23, 2);
		DFTOutMem[i] = val;

	}
	fft_config2.setDir(false, 0);
	hls::fft<config2_recon>(DFTOutMem, IMOutMem, &fft_status2, &fft_config2);
	//ifft_out_t sigouttemp[MEMSIZE];

	for(int i=0;i<MEMSIZE;i++) {
#pragma HLS PIPELINE
		cmpxIFFTOut datOut = IMOutMem[i];
		ifft_out_t  datV = datOut.real();
		//sigOut.write( datV);
		sigOut[i] = datV;
	}

}

#define HLS_8U_MIN   0
#define HLS_8U_MAX   255
template<typename T>
unsigned char Convert2uchar(T v)
{
    unsigned char result=HLS_8U_MIN;
    if(v>=HLS_8U_MAX)
    {
        result=HLS_8U_MAX;
    }
    else if(v>=HLS_8U_MIN&&v<HLS_8U_MAX)
    {
        ap_fixed<9,9,AP_RND> temp=v;
        result=temp;
    }
    else if(v<HLS_8U_MIN){
    	ap_fixed<9,9,AP_RND> temp = v*-1;
    	result = temp;
    }
    return result;
}
void pyr_recon_ifft(hls::stream< t_recon_complex>  &DFTOut,  hls::stream<  PIXEL_RAW>  &RGBOut,hls::stream<ifft_out_t> &outf) {

#pragma HLS data_pack variable=DFTOut
#pragma HLS data_pack variable=RGBOut

	//hls::stream< ifft_out_t>  sigOut;
	ifft_out_t sigOut[MEMSIZE];
	t_recon_complex DFTin[MEMSIZE];
#pragma HLS STREAM variable=DFTin depth=512 dim=1
#pragma HLS STREAM variable=sigOut depth=512 dim=1
#pragma HLS DATAFLOW
	for(int i=0;i<MEMSIZE;i++) {
#pragma HLS PIPELINE
		DFTin[i] = DFTOut.read();
	}
	cmifft(DFTin, sigOut);
	//cmifft(DFTOut[1], sigOut[1]);
	//cmifft(DFTOut[2], sigOut[2]);

	for(int i=0;i<512;i++ ){
#pragma HLS INLINE
		ifft_out_t tRs = sigOut[i]+2 ;
		ifft_out_t Rs = tRs<<6;
		outf.write(tRs);

		//ifft_out_t Gs = sigOut[1].read() <<10;
		//ifft_out_t Bs = sigOut[2].read() <<10;
		//ap_fixed<
		unsigned char R;// ,G,B;
		R = Convert2uchar(Rs);
		//G = Convert2uchar(Gs);
		//B = Convert2uchar(Bs);
		//printf("%f  %f  %d\n",sigOut[i].to_float(),Rs.to_float(),R);
		PIXEL_RAW p;
		p.range(7,0) = R;
		p.range(15,8) = 0;
		p.range(23,16)  = 0;

		RGBOut.write(p );
	}

}

//
//void pyr_recon_combine(hls::stream< cmpxDataOut>  &pyrFFT,
//				  hls::stream< PIXEL_RAW>   &sigOut) {
//
//	#pragma HLS data_pack variable=pyrFFT
//	#pragma HLS data_pack variable=sigOut
//
//
//
//	#pragma HLS DATAFLOW
//	hls::stream< t_recon_complex>  DFTOut[reconC];
//	pyr_recon(pyrFFT, DFTOut);
//	pyr_recon_ifft(DFTOut, sigOut);
//
//}


