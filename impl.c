#ifndef _MUTIPLY_H
#define _MULTIPLY_H

void naive_multiply(int *src1, int *src2, int *dst, int src1_w, int src1_h,
                    int src2_w, int src2_h)
{
    for (int i = 0; i < src1_h; ++i) {
        for (int j = 0; j < src2_w; ++j) {
            dst[i * src2_w + j] = 0;
            for (int k = 0; k < src2_h; ++k) {
                dst[i * src2_w + j] += src1[i * src1_w + k] * src2[k * src2_w + j];
            }
        }
    }
}

void sse_multiply(int *src1, int *src2, int *dst, int src1_w, int src1_h,
                  int src2_w, int src2_h)
{
    for (int x = 0; x < src2_h; x += 4) {
        for (int y = 0; y < src2_w; y += 4) {
        
            __m128i I0 = _mm_loadu_si128((__m128i *)(src2 + (y + 0) * src2_w + x));
            __m128i I1 = _mm_loadu_si128((__m128i *)(src2 + (y + 1) * src2_w + x));
            __m128i I2 = _mm_loadu_si128((__m128i *)(src2 + (y + 2) * src2_w + x));
            __m128i I3 = _mm_loadu_si128((__m128i *)(src2 + (y + 3) * src2_w + x));

            __m128i I4 = _mm_loadu_si128((__m128i *)(src1 +0*src1_w));
            __m128i I5 = _mm_loadu_si128((__m128i *)(src1 +1*src1_w));
            __m128i I6 = _mm_loadu_si128((__m128i *)(src1 +2*src1_w));
            __m128i I7 = _mm_loadu_si128((__m128i *)(src1 +3*src1_w));  
           
            __m128i T0 = _mm_unpacklo_epi32(I0, I1);
            __m128i T1 = _mm_unpacklo_epi32(I2, I3);
            __m128i T2 = _mm_unpackhi_epi32(I0, I1);
            __m128i T3 = _mm_unpackhi_epi32(I2, I3);
           
            I0 = _mm_unpacklo_epi64(T0, T1);
            I1 = _mm_unpackhi_epi64(T0, T1);
            I2 = _mm_unpacklo_epi64(T2, T3);
            I3 = _mm_unpackhi_epi64(T2, T3);
/****** the first row ******/	   
	    __m128i T40 =  _mm_mullo_epi16(I4,I0);  
            __m128i T41 =  _mm_mullo_epi16(I4,I1);  
            __m128i T42 =  _mm_mullo_epi16(I4,I2);  
	    __m128i T43 =  _mm_mullo_epi16(I4,I3);  

            __m128i T00 = _mm_unpacklo_epi32(T40, T41);
            __m128i T10 = _mm_unpacklo_epi32(T42, T43);
            __m128i T20 = _mm_unpackhi_epi32(T40, T41);
            __m128i T30 = _mm_unpackhi_epi32(T42, T43);
           
            T40 = _mm_unpacklo_epi64(T00, T10);
            T41 = _mm_unpackhi_epi64(T00, T10);
            T42 = _mm_unpacklo_epi64(T20, T30);
            T43 = _mm_unpackhi_epi64(T20, T30);
 
            __m128i A0 = _mm_add_epi64(T40, T41);
            __m128i A1 = _mm_add_epi64(T42, T43);
            __m128i AA = _mm_add_epi64(A0, A1);
/****** the second row ******/	
            __m128i T50 =  _mm_mullo_epi16(I5,I0);  
            __m128i T51 =  _mm_mullo_epi16(I5,I1);  
            __m128i T52 =  _mm_mullo_epi16(I5,I2);  
            __m128i T53 =  _mm_mullo_epi16(I5,I3);

            __m128i T01 = _mm_unpacklo_epi32(T50, T51);
            __m128i T11 = _mm_unpacklo_epi32(T52, T53);
            __m128i T21 = _mm_unpackhi_epi32(T50, T51);
            __m128i T31 = _mm_unpackhi_epi32(T52, T53);
           
            T50 = _mm_unpacklo_epi64(T01, T11);
            T51 = _mm_unpackhi_epi64(T01, T11);
            T52 = _mm_unpacklo_epi64(T21, T31);
            T53 = _mm_unpackhi_epi64(T21, T31);

            __m128i B0 = _mm_add_epi64(T50, T51);
            __m128i B1 = _mm_add_epi64(T52, T53);
            __m128i BB = _mm_add_epi64(B0, B1);
/****** the third row *****/ 
	    __m128i T60 =  _mm_mullo_epi16(I6,I0);
            __m128i T61 =  _mm_mullo_epi16(I6,I1);
            __m128i T62 =  _mm_mullo_epi16(I6,I2);
            __m128i T63 =  _mm_mullo_epi16(I6,I3);

            __m128i T02 = _mm_unpacklo_epi32(T60, T61);
            __m128i T12 = _mm_unpacklo_epi32(T62, T63);
            __m128i T22 = _mm_unpackhi_epi32(T60, T61);
            __m128i T32 = _mm_unpackhi_epi32(T62, T63);
           
            T60 = _mm_unpacklo_epi64(T02, T12);
            T61 = _mm_unpackhi_epi64(T02, T12);
            T62 = _mm_unpacklo_epi64(T22, T32);
            T63 = _mm_unpackhi_epi64(T22, T32);

            __m128i C0 = _mm_add_epi64(T60, T61);
            __m128i C1 = _mm_add_epi64(T62, T63);
            __m128i CC = _mm_add_epi64(C0, C1);  
/****** the forth row ******/
            __m128i T70 =  _mm_mullo_epi16(I7,I0);
            __m128i T71 =  _mm_mullo_epi16(I7,I1);
            __m128i T72 =  _mm_mullo_epi16(I7,I2);
            __m128i T73 =  _mm_mullo_epi16(I7,I3);

	    __m128i T03 = _mm_unpacklo_epi32(T70, T71);
            __m128i T13 = _mm_unpacklo_epi32(T72, T73);
            __m128i T23 = _mm_unpackhi_epi32(T70, T71);
            __m128i T33 = _mm_unpackhi_epi32(T72, T73);
           
            T70 = _mm_unpacklo_epi64(T03, T13);
            T71 = _mm_unpackhi_epi64(T03, T13);
            T72 = _mm_unpacklo_epi64(T23, T33);
            T73 = _mm_unpackhi_epi64(T23, T33);

            __m128i D0 = _mm_add_epi64(T70, T71);
            __m128i D1 = _mm_add_epi64(T72, T73);
            __m128i DD = _mm_add_epi64(D0, D1);
           
	    _mm_storeu_si128((__m128i *)(dst + ((x + 0) * 4) + y), AA);
            _mm_storeu_si128((__m128i *)(dst + ((x + 1) * 4) + y), BB);
            _mm_storeu_si128((__m128i *)(dst + ((x + 2) * 4) + y), CC);
            _mm_storeu_si128((__m128i *)(dst + ((x + 3) * 4) + y), DD);
    
	}
    }	
      

}

void sse_prefetch_multiply(int *src1, int *src2, int *dst, int src1_w,
                           int src1_h, int src2_w, int src2_h)
{
}

#endif
