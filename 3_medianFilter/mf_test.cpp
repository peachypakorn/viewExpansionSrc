#include "median_filter.h"
#include <stdio.h>
#include <string.h>

int main() {
int width = 1520;
int hight =10;
hls::stream<data_t> output("test");
hls::stream<data_t> src("test2");

int i ;

for(int j = 0;j<10 ;j++){
for (i = 0; i < 1520; i++) {

	data_t val;
	val = (data_t)temp[i + j*1520];
	//else if(j==1)val = (t_output_scalar)temp2[i];
	//else if(j==2)val = (t_output_scalar)temp3[i];

	//val.last =0;
	src.write(val);
	//printf("%d  %f \n",i,temp[i]);
}
}
median_strm(width,hight,src,output);

FILE * fo = fopen("mdf_out.txt", " wb");
for (int j = 0; j <10; ++j) {

for( i=0;i<1520;i++)
	{

		data_t val2 = output.read();
		//data_t val = 0.789;
		//std::cout << val << std::endl;
		fprintf(fo, "%f  , ", val2.to_double());
		//printf("%f \n",val.data.to_double());
		//printf("%d \n",i);
	}
fprintf(fo," \n");
}

	fclose(fo);
	return 0;
}

