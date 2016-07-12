#include <assert.h>
#include <stdint.h>
#include <hls_stream.h>
#include "hls_dsp.h"
#include "median_filter.h"

template <typename T, int KMED>
T median(T window[KMED * KMED]) {
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_RESHAPE variable=window complete dim=1

	int const N=KMED*KMED;
	T t[KMED * KMED], z[KMED * KMED];


//#pragma HLS ARRAY_RESHAPE variable=t complete dim=1
//#pragma HLS ARRAY_RESHAPE variable=z complete dim=1

#pragma HLS	ARRAY_PARTITION variable=t complete dim=0
#pragma HLS	ARRAY_PARTITION variable=z complete dim=0

	unsigned char i, k, stage;
	// Copy input
	for(i=0;i<N;i++)  z[i] = window[i];
	// Sorting network loop
	for(stage = 1; stage <= N; stage++ )
		{
			//if ((stage%2)==1) k = 0;
			//if ((stage%2)==0) k = 1;
			k = (stage%2==1) ? 0 : 1;
			for(i=k;i<N-1;i=i+2) {
				float a = z[i];
				float b = z[i+1];
				t[i] = MIN(z[i], z[i+1]);
				float min =t[i];
				t[i+1] = MAX(z[i], z[i+1]);
				float max = t[i+1];
				z[i] = t[i];
				z[i+1] = t[i+1];
			}
		}
float finish =z[N/2];
	return z[N/2];//return Medien
}


template <typename T, int K>
void median_str(int width, int height,
        hls::stream<T> &src, hls::stream<T> &dst) {
    T window[ K * K], pix, pixel[K], med;
    const int K2 = K/2;

	// Line-buffers allowing full pixel reuse in vertical pass
    static T line_buffer[K][MAX_WIDTH];
#pragma HLS DATAFLOW
#pragma HLS ARRAY_PARTITION variable=line_buffer dim=1 complete
#pragma HLS INLINE // Into a DATAFLOW region
    // These assertions let HLS know the upper bounds of loops
    assert(height < MAX_IMG_ROWS); //1080
    assert(width < MAX_IMG_COLS);  //1920
//    assert(vconv_xlim < MAX_IMG_COLS - (K - 1));
    short int r,c;

L1:for(r=0; r<height; r++) {
#pragma HLS LOOP_TRIPCOUNT min=600 max=1080 avg=720
    	L2:for(c=0;c<width;c++) {
#pragma HLS LOOP_TRIPCOUNT min=512 max=1920 avg=1600
#pragma HLS PIPELINE II=1
    		for(int i=0;i<K-1;i++) {
    			line_buffer[i][c] = line_buffer[i+1][c];
    			pixel[i] = line_buffer[i][c];
    		}
    		pix = src.read();
    		pixel[K -1] = line_buffer[K-1][c]=pix;

    		for(int i=0;i<K;i++)
    			for(int j=0;j<K-1;j++) {
    				window[i*K +j] = window[i*K +j+1];
    			}
    		for(int i=0;i<K;i++)
    			window[i*K + K - 1] = pixel[i];
    		if (r>=K && c>=K) {
    			med = median<T, K>(window);
    			pix = med;
    		}
    		dst.write(pix);
    		//if (r>0 && c>0)

    	}
    }

}


void median_strm(
        int width, int height,
        hls::stream<data_t> &src, hls::stream<data_t> &dst) {
//#pragma HLS INTERFACE ap_none port=height metadata="-bus_bundle hls_ctrl"
//#pragma HLS INTERFACE ap_none port=width metadata="-bus_bundle hls_ctrl"
//#pragma HLS INTERFACE s_axilite port=return metadata="-bus_bundle hls_ctrl"
#pragma HLS INTERFACE port=src axis depth=6080 // TEST_IMG_SIZE
#pragma HLS INTERFACE port=dst axis depth=6080 // TEST_IMG_SIZE

//#pragma HLS RESOURCE variable=width      core=AXI4LiteS metadata="-bus_bundle hls_ctrl"
//#pragma HLS RESOURCE variable=height      core=AXI4LiteS metadata="-bus_bundle hls_ctrl"

//#pragma HLS RESOURCE variable=return core=AXI4LiteS metadata="-bus_bundle hls_ctrl"
//#pragma HLS RESOURCE variable=width      core=AXI4LiteS metadata="-bus_bundle hls_ctrl"
//#pragma HLS RESOURCE variable=height      core=AXI4LiteS metadata="-bus_bundle hls_ctrl"

//#pragma HLS DATAFLOW
//#pragma HLS INLINE region // bring loops in sub-functions to this DATAFLOW region

	median_str<data_t, 3>(width, height,     src, dst);

}
