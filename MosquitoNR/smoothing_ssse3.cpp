//------------------------------------------------------------------------------
//		smoothing_ssse3.cpp
//------------------------------------------------------------------------------

#include "mosquito_nr.h"

#if !defined(_WIN64)
#define rax	eax
#define rbx	ebx
#define rcx	ecx
#define rdx	edx
#define rsi	esi
#define rdi	edi
#define rbp	ebp
#else
#define rax	rax
#define rbx	rbx
#define rcx	rcx
#define rdx	rdx
#define rsi	rsi
#define rdi	rdi
#define rbp	rbp
#endif

// direction-aware blur
void MosquitoNR::SmoothingSSSE3(int thread_id)
{
	const int y_start = height *  thread_id      / threads;
	const int y_end   = height * (thread_id + 1) / threads;
	if (y_start == y_end) return;
	const int width  = this->width;
	const int pitch  = this->pitch;
	const int pitch2 = pitch * sizeof(short);
	__declspec(align(16)) short sad[48], tmp[16];
	short *srcp, *dstp, *sadp = sad, *tmpp = tmp;
	for (int i = 0; i <  8; ++i) tmp[i] = 4;
	for (int i = 8; i < 16; ++i) tmp[i] = 3;

	if (radius == 1)
	{
		const int coef0 =  64 - strength * 2;	// own pixel's coefficient (when divisor = 64)
		const int coef1 = 128 - strength * 4;	// own pixel's coefficient (when divisor = 128)
		const int coef2 = strength;				// other pixel's coefficient

		for (int y = y_start; y < y_end; ++y)
		{
			srcp = luma[0] + (y + 2) * pitch + 8;
			dstp = luma[1] + (y + 2) * pitch + 8;

			for (int x = 0; x < width; x += 8)
			{
				__asm
				{
					mov			rsi, srcp
					mov			rdi, sadp				// edi = sadp
					mov			rdx, tmpp				// edx = tmpp
					mov			eax, pitch2				// eax = pitch * sizeof(short)
					movdqa		xmm6, [rdx]				// xmm6 = [4] * 8
					movdqa		xmm7, [rsi]				// xmm7 = (  0,  0 )
					sub			rsi, rax				// esi = srcp - pitch

					movdqu		xmm0, [rsi+rax-2]		// xmm0 = ( -1,  0 )
					movdqu		xmm1, [rsi+rax+2]		// xmm1 = (  1,  0 )
					movdqu		xmm2, [rsi-2]			// xmm2 = ( -1, -1 )
					movdqu		xmm3, [rsi+2*rax+2]		// xmm3 = (  1,  1 )

					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					paddw		xmm0, xmm1
					movdqa		[rdi], xmm0

					movdqa		xmm0, [rsi]				// xmm0 = (  0, -1 )
					movdqa		xmm1, [rsi+2*rax]		// xmm1 = (  0,  1 )

					paddw		xmm4, xmm2
					paddw		xmm5, xmm3
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm4, xmm5
					paddw		xmm4, xmm6				// add "identification number" to the lower 3 bits (4)
					psubw		xmm6, [rdx+16]			// (The lower 3 bits are always zero.
					movdqa		[rdi+16], xmm4			//  If input is 9-bit or more, this hack doesn't work)

					movdqa		xmm4, xmm2
					movdqa		xmm5, xmm3
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm2, xmm3
					paddw		xmm2, xmm6				// (1)
					paddw		xmm6, [rdx]
					movdqa		[rdi+32], xmm2

					movdqu		xmm2, [rsi+2]			// xmm2 = (  1, -1 )
					movdqu		xmm3, [rsi+2*rax-2]		// xmm3 = ( -1,  1 )

					paddw		xmm4, xmm0
					paddw		xmm5, xmm1
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm4, xmm5
					paddw		xmm4, xmm6				// (5)
					psubw		xmm6, [rdx+16]
					movdqa		[rdi+48], xmm4

					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					paddw		xmm0, xmm1
					paddw		xmm0, xmm6				// (2)
					paddw		xmm6, [rdx]
					movdqa		[rdi+64], xmm0

					movdqu		xmm0, [rsi+rax+2]		// xmm0 = (  1,  0 )
					movdqu		xmm1, [rsi+rax-2]		// xmm1 = ( -1,  0 )

					paddw		xmm4, xmm2
					paddw		xmm5, xmm3
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm4, xmm5
					paddw		xmm4, xmm6				// (6)
					psubw		xmm6, [rdx+16]

					pminsw		xmm4, [rdi]
					pminsw		xmm4, [rdi+16]
					pminsw		xmm4, [rdi+32]
					pminsw		xmm4, [rdi+48]
					pminsw		xmm4, [rdi+64]

					paddw		xmm0, xmm2
					paddw		xmm1, xmm3
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm2, xmm3
					paddw		xmm2, xmm6				// (3)
					paddw		xmm6, [rdx]

					psraw		xmm0, 1
					psraw		xmm1, 1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					paddw		xmm0, xmm1
					paddw		xmm0, xmm6				// (7)

					pminsw		xmm4, xmm2
					pminsw		xmm4, xmm0
					movdqa		[rdi], xmm4
				}
				
				for (int i = 0; i < 8; ++i, ++srcp, ++dstp)
				{
					if ((sad[i] &~ 7) == 0) { *dstp = *srcp; continue; }

					switch (sad[i] & 7)
					{
						case 0:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-1]       + srcp[1]      ) + 32) >> 6; break;
						case 1:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch-1] + srcp[pitch+1]) + 32) >> 6; break;
						case 2:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch]   + srcp[pitch]  ) + 32) >> 6; break;
						case 3:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch+1] + srcp[pitch-1]) + 32) >> 6; break;
						case 4:
							*dstp = (coef1 * srcp[0] + coef2 * (srcp[-pitch-1] + srcp[-1]     + srcp[1]     + srcp[pitch+1]) + 64) >> 7; break;
						case 5:
							*dstp = (coef1 * srcp[0] + coef2 * (srcp[-pitch-1] + srcp[-pitch] + srcp[pitch] + srcp[pitch+1]) + 64) >> 7; break;
						case 6:
							*dstp = (coef1 * srcp[0] + coef2 * (srcp[-pitch+1] + srcp[-pitch] + srcp[pitch] + srcp[pitch-1]) + 64) >> 7; break;
						case 7:
							*dstp = (coef1 * srcp[0] + coef2 * (srcp[-pitch+1] + srcp[1]      + srcp[-1]    + srcp[pitch-1]) + 64) >> 7; break;
					}
				}
			}
		}
	}
	else	// radius == 2
	{
		const int coef0 = 128 - strength * 4;	// own pixel's coefficient (when divisor = 128)
		const int coef1 = 256 - strength * 8;	// own pixel's coefficient (when divisor = 256)
		const int coef2 = strength;				// other pixel's coefficient
		const int coef3 = strength * 2;			// other pixel's coefficient (doubled)

		for (int y = y_start; y < y_end; ++y)
		{
			srcp = luma[0] + (y + 2) * pitch + 8;
			dstp = luma[1] + (y + 2) * pitch + 8;

			for (int x = 0; x < width; x += 8)
			{
				__asm
				{
					mov			rsi, srcp
					mov			rdi, sadp				// edi = sadp
					mov			rdx, tmpp				// edx = tmpp
					mov			eax, pitch2				// eax = pitch * sizeof(short)
					lea			ecx, [eax+2*eax]		// ecx = pitch * sizeof(short) * 3
					movdqa		xmm6, [rdx]				// xmm6 = [0x0001] * 8
					movdqa		xmm7, [rsi]				// xmm7 = (  0,  0 )
					sub			rsi, rax
					sub			rsi, rax				// esi = srcp - 2 * pitch

					movdqu		xmm0, [rsi+2*rax-2]		// xmm0 = ( -1,  0 )
					movdqu		xmm1, [rsi+2*rax+2]		// xmm1 = (  1,  0 )
					movdqu		xmm2, [rsi+2*rax-4]		// xmm2 = ( -2,  0 )
					movdqu		xmm3, [rsi+2*rax+4]		// xmm3 = (  2,  0 )
					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm0, xmm1
					paddw		xmm2, xmm3
					paddw		xmm0, xmm2
					movdqa		[rdi], xmm0

					movdqu		xmm0, [rsi+rax-2]		// xmm0 = ( -1, -1 )
					movdqu		xmm1, [rsi+rcx+2]		// xmm1 = (  1,  1 )
					movdqu		xmm2, [rsi+rax-4]		// xmm2 = ( -2, -1 )
					movdqu		xmm3, [rsi+rcx+4]		// xmm3 = (  2,  1 )
					paddw		xmm4, xmm0
					paddw		xmm5, xmm1
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm2, xmm3
					paddw		xmm4, xmm5
					paddw		xmm2, xmm4
					paddw		xmm2, xmm6				// add "identification number" to the lower 3 bits (4)
					psubw		xmm6, [rdx+16]
					movdqa		[rdi+16], xmm2

					movdqu		xmm2, [rsi-4]			// xmm2 = ( -2, -2 )
					movdqu		xmm3, [rsi+4*rax+4]		// xmm3 = (  2,  2 )
					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm0, xmm1
					paddw		xmm2, xmm3
					paddw		xmm0, xmm2
					paddw		xmm0, xmm6				// (1)
					paddw		xmm6, [rdx]
					movdqa		[rdi+32], xmm0

					movdqa		xmm0, [rsi+rax]			// xmm0 = (  0, -1 )
					movdqa		xmm1, [rsi+rcx]			// xmm1 = (  0,  1 )
					movdqu		xmm2, [rsi-2]			// xmm2 = ( -1, -2 )
					movdqu		xmm3, [rsi+4*rax+2]		// xmm3 = (  1,  2 )
					paddw		xmm4, xmm0
					paddw		xmm5, xmm1
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm2, xmm3
					paddw		xmm4, xmm5
					paddw		xmm2, xmm4
					paddw		xmm2, xmm6				// (5)
					psubw		xmm6, [rdx+16]
					movdqa		[rdi+48], xmm2

					movdqa		xmm2, [rsi]				// xmm2 = (  0, -2 )
					movdqa		xmm3, [rsi+4*rax]		// xmm3 = (  0,  2 )
					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm0, xmm1
					paddw		xmm2, xmm3
					paddw		xmm0, xmm2
					paddw		xmm0, xmm6				// (2)
					paddw		xmm6, [rdx]
					movdqa		[rdi+64], xmm0

					movdqu		xmm0, [rsi+rax+2]		// xmm0 = (  1, -1 )
					movdqu		xmm1, [rsi+rcx-2]		// xmm1 = ( -1,  1 )
					movdqu		xmm2, [rsi+2]			// xmm2 = (  1, -2 )
					movdqu		xmm3, [rsi+4*rax-2]		// xmm3 = ( -1,  2 )
					paddw		xmm4, xmm0
					paddw		xmm5, xmm1
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm2, xmm3
					paddw		xmm4, xmm5
					paddw		xmm2, xmm4
					paddw		xmm2, xmm6				// (6)
					psubw		xmm6, [rdx+16]
					movdqa		[rdi+80], xmm2

					movdqu		xmm2, [rsi+4]			// xmm2 = (  2, -2 )
					movdqu		xmm3, [rsi+4*rax-4]		// xmm3 = ( -2,  2 )
					movdqa		xmm4, xmm0
					movdqa		xmm5, xmm1
					psubw		xmm0, xmm7
					psubw		xmm1, xmm7
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					pabsw		xmm0, xmm0
					pabsw		xmm1, xmm1
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					paddw		xmm0, xmm1
					paddw		xmm2, xmm3
					paddw		xmm0, xmm2
					paddw		xmm0, xmm6				// (3)
					paddw		xmm6, [rdx]

					pminsw		xmm0, [rdi]
					pminsw		xmm0, [rdi+16]
					pminsw		xmm0, [rdi+32]
					pminsw		xmm0, [rdi+48]
					pminsw		xmm0, [rdi+64]
					pminsw		xmm0, [rdi+80]

					movdqu		xmm1, [rsi+2*rax+2]		// xmm1 = (  1,  0 )
					movdqu		xmm2, [rsi+rax+4]		// xmm2 = (  2, -1 )
					movdqu		xmm3, [rsi+rcx-4]		// xmm3 = ( -2,  1 )
					paddw		xmm4, xmm1
					movdqu		xmm1, [rsi+2*rax-2]		// xmm1 = ( -1,  0 )
					paddw		xmm5, xmm1
					psraw		xmm4, 1
					psraw		xmm5, 1
					psubw		xmm2, xmm7
					psubw		xmm3, xmm7
					psubw		xmm4, xmm7
					psubw		xmm5, xmm7
					pabsw		xmm2, xmm2
					pabsw		xmm3, xmm3
					pabsw		xmm4, xmm4
					pabsw		xmm5, xmm5
					paddw		xmm2, xmm3
					paddw		xmm4, xmm5
					paddw		xmm2, xmm4
					paddw		xmm2, xmm6				// (7)

					pminsw		xmm0, xmm2
					movdqa		[rdi], xmm0
				}
			
				for (int i = 0; i < 8; ++i, ++srcp, ++dstp)
				{
					if ((sad[i] &~ 7) == 0) { *dstp = *srcp; continue; }

					switch (sad[i] & 7)
					{
						case 0:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-2]        + srcp[-1]       + srcp[1]       + srcp[2]       ) + 64) >> 7; break;
						case 1:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch2-2] + srcp[-pitch-1] + srcp[pitch+1] + srcp[pitch2+2]) + 64) >> 7; break;
						case 2:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch2]   + srcp[-pitch]   + srcp[pitch]   + srcp[pitch2]  ) + 64) >> 7; break;
						case 3:
							*dstp = (coef0 * srcp[0] + coef2 * (srcp[-pitch2+2] + srcp[-pitch+1] + srcp[pitch-1] + srcp[pitch2-2]) + 64) >> 7; break;
						case 4:
							*dstp = (coef1 * srcp[0] + coef3 * (srcp[-pitch -2] + srcp[pitch +2]) + coef2 * (srcp[-pitch-1] + srcp[-1]     + srcp[1]     + srcp[pitch+1]) + 128) >> 8; break;
						case 5:
							*dstp = (coef1 * srcp[0] + coef3 * (srcp[-pitch2-1] + srcp[pitch2+1]) + coef2 * (srcp[-pitch-1] + srcp[-pitch] + srcp[pitch] + srcp[pitch+1]) + 128) >> 8; break;
						case 6:
							*dstp = (coef1 * srcp[0] + coef3 * (srcp[-pitch2+1] + srcp[pitch2-1]) + coef2 * (srcp[-pitch+1] + srcp[-pitch] + srcp[pitch] + srcp[pitch-1]) + 128) >> 8; break;
						case 7:
							*dstp = (coef1 * srcp[0] + coef3 * (srcp[-pitch +2] + srcp[pitch -2]) + coef2 * (srcp[-pitch+1] + srcp[1]      + srcp[-1]    + srcp[pitch-1]) + 128) >> 8; break;
					}
				}
			}
		}
	}

	// vertical reflection
	if (y_start <= 1 && 1 < y_end)
		memcpy(luma[1] + pitch, luma[1] + 3 * pitch, pitch * sizeof(short));
	if (y_start <= 2 && 2 < y_end)
		memcpy(luma[1],         luma[1] + 4 * pitch, pitch * sizeof(short));
	if (y_start <= height - 3 && height - 3 < y_end)
		memcpy(luma[1] + (height + 3) * pitch, luma[1] + (height - 1) * pitch, pitch * sizeof(short));
	if (y_start <= height - 2 && height - 2 < y_end)
		memcpy(luma[1] + (height + 2) * pitch, luma[1] +  height      * pitch, pitch * sizeof(short));
}
