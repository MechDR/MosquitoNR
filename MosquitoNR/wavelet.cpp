//------------------------------------------------------------------------------
//		wavelet.cpp
//------------------------------------------------------------------------------

/*
	To divide an image into low and high frequency components, CDF 5/3 wavelet is used.
	Wavelet transform is applied to rows and columns independently.
	Outside the image, reflected data must be prepared. (like 0123.. -> 210123..)

	<Forward transform>
	1. From odd number pixels, subtract an average of both neighboring pixel.
	   This data is called detail coefficients, and shows high frequency components.
	   ex) [1] -= ([0] + [2]) / 2
	2. To even number pixels, add a quarter of a sum of both neighboring pixel (i.e. detail coefficients).
	   This data is called approximation coefficients, and shows low frequency components.
	   ex) [0] += ([1] + [1]) / 4
	       [2] += ([1] + [3]) / 4

	<Inverse transform>
	Just do opposite operations to the forward transform.

	CDF 5/3 wavelet transform is reversible including edge points.
*/

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

void MosquitoNR::WaveletVert1(int thread_id)
{
	const int y_start = (height + 7) / 8 *  thread_id      / threads * 8;
	const int y_end   = (height + 7) / 8 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width = this->width;
	const int pitch = this->pitch;
	const int hloop = (width + 7) / 8;

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp = luma[0] + y * pitch + 8;
		short* dstp = bufy[0] + y / 2 * pitch + 8;

		__asm
		{
			mov			rsi, srcp			// esi = srcp
			mov			rdi, dstp			// edi = dstp
			mov			eax, pitch
			mov			ecx, hloop			// ecx = hloop
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3

align 16
next8columns:
#if !defined(_WIN64)
			push		esi
#else
			mov			r12, rsi
#endif
			movdqa		xmm2, [rsi]
			movdqa		xmm0, [rsi+rax]
			movdqa		xmm1, [rsi+2*rax]
			paddw		xmm2, xmm1
			psraw		xmm2, 1
			psubw		xmm0, xmm2
			add			esi, ebx

			movdqa		xmm2, [rsi]
			movdqa		xmm3, [rsi+rax]
			movdqa		xmm4, [rsi+2*rax]
			movdqa		xmm5, [rsi+rbx]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 2
			psraw		xmm2, 2
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi], xmm6
			movdqa		[rdi+rax], xmm7
			lea			rsi, [rsi+4*rax]

			movdqa		xmm0, [rsi]
			movdqa		xmm1, [rsi+rax]
			movdqa		xmm2, [rsi+2*rax]
			movdqa		xmm3, [rsi+rbx]
			movdqa		xmm6, xmm5
			movdqa		xmm7, xmm1
			paddw		xmm5, xmm1
			paddw		xmm1, xmm3
			psraw		xmm5, 1
			psraw		xmm1, 1
			psubw		xmm0, xmm5
			psubw		xmm2, xmm1
			paddw		xmm4, xmm0
			paddw		xmm0, xmm2
			psraw		xmm4, 2
			psraw		xmm0, 2
			paddw		xmm6, xmm4
			paddw		xmm7, xmm0
			movdqa		[rdi+2*rax], xmm6
			movdqa		[rdi+rbx], xmm7

#if !defined(_WIN64)
			pop			esi
#else
			mov			rsi, r12
#endif
			add			rsi, 16
			add			rdi, 16
			sub			ecx, 1
			jnz			next8columns
		}

		// horizontal reflection
		short* p = dstp;
		for (int i = 0; i < 4; ++i, p += pitch)
			p[-2] = p[2], p[-1] = p[1], p[width] = p[width-2], p[width+1] = p[width-3];
	}
}

void MosquitoNR::WaveletHorz1(int thread_id)
{
	const int y_start = (height + 15) / 16 *  thread_id      / threads * 8;
	const int y_end   = (height + 15) / 16 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width  = this->width;
	const int pitch  = this->pitch;
	const int hloop1 = (width + 4 + 2 + 3) / 4;
	const int hloop2 = (width + 3) / 4;
	short* work = this->work[thread_id];

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp = bufy[0] + y * pitch + 4;
		short* dstp = luma[0] + y / 2 * pitch + 8;

		__asm
		{
			// shuffle
			mov			rsi, srcp			// esi = srcp
			mov			rdi, work			// edi = work
			mov			eax, pitch
			mov			ecx, hloop1			// ecx = hloop1
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3
			lea			rdx, [rsi+4*rax]	// edx = srcp + 4 * pitch

align 16
shuffle_next4columns:
			movq		xmm0, qword ptr [rsi]				// 03, 02, 01, 00
			movq		xmm1, qword ptr [rsi+rax]			// 13, 12, 11, 10
			movq		xmm2, qword ptr [rsi+2*rax]			// 23, 22, 21, 20
			movq		xmm3, qword ptr [rsi+rbx]			// 33, 32, 31, 30
			movq		xmm4, qword ptr [rdx]				// 43, 42, 41, 40
			movq		xmm5, qword ptr [rdx+rax]			// 53, 52, 51, 50
			movq		xmm6, qword ptr [rdx+2*rax]			// 63, 62, 61, 60
			movq		xmm7, qword ptr [rdx+rbx]			// 73, 72, 71, 70
			punpcklwd	xmm0, xmm1			// 13, 03, 12, 02, 11, 01, 10, 00
			punpcklwd	xmm2, xmm3			// 33, 23, 32, 22, 31, 21, 30, 20
			punpcklwd	xmm4, xmm5			// 53, 43, 52, 42, 51, 41, 50, 40
			punpcklwd	xmm6, xmm7			// 73, 63, 72, 62, 71, 61, 70, 60
			movdqa		xmm1, xmm0
			movdqa		xmm5, xmm4
			punpckldq	xmm0, xmm2			// 31, 21, 11, 01, 30, 20, 10, 00
			punpckhdq	xmm1, xmm2			// 33, 23, 13, 03, 32, 22, 12, 02
			punpckldq	xmm4, xmm6			// 71, 61, 51, 41, 70, 60, 50, 40
			punpckhdq	xmm5, xmm6			// 73, 63, 53, 43, 72, 62, 52, 42
			movdqa		xmm2, xmm0
			movdqa		xmm3, xmm1
			punpcklqdq	xmm0, xmm4			// 70, 60, 50, 40, 30, 20, 10, 00
			punpckhqdq	xmm2, xmm4			// 71, 61, 51, 41, 31, 21, 11, 01
			punpcklqdq	xmm1, xmm5			// 72, 62, 52, 42, 32, 22, 12, 02
			punpckhqdq	xmm3, xmm5			// 73, 63, 53, 43, 33, 23, 13, 03
			movdqa		[rdi], xmm0
			movdqa		[rdi+16], xmm2
			movdqa		[rdi+32], xmm1
			movdqa		[rdi+48], xmm3
			add			rsi, 8
			add			rdx, 8
			add			rdi, 64
			sub			ecx, 1
			jnz			shuffle_next4columns

			// wavelet transform
			mov			rsi, work
			mov			rdi, dstp			// edi = dstp
			mov			ecx, hloop2			// ecx = hloop2
			add			rsi, 64				// esi = work + 32

			movdqa		xmm2, [rsi-32]
			movdqa		xmm0, [rsi-16]
			movdqa		xmm1, [rsi]
			paddw		xmm2, xmm1
			psraw		xmm2, 1
			psubw		xmm0, xmm2

align 16
wavelet_next4columns:
			movdqa		xmm2, [rsi+16]
			movdqa		xmm3, [rsi+32]
			movdqa		xmm4, [rsi+48]
			movdqa		xmm5, [rsi+64]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 2
			psraw		xmm2, 2
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi], xmm6
			movdqa		[rdi+16], xmm7
			movdqa		xmm0, xmm4
			movdqa		xmm1, xmm5
			add			rsi, 64
			add			rdi, 32
			add			rdx, 32
			sub			ecx, 1
			jnz			wavelet_next4columns

			// horizontal reflection
			mov			eax, width
			test		eax, 1
			jnz			no_reflection		// if (width % 2 == 0) {
			mov			rdi, dstp			//   edi = dstp
			shl			eax, 3				//   eax = width / 2 * 16
			movdqa		xmm0, [rdi+rax-16]
			movdqa		[rdi+rax], xmm0
no_reflection:								// }
		}
	}
}

void MosquitoNR::WaveletVert2(int thread_id)
{
	const int y_start = (height + 7) / 8 *  thread_id      / threads * 8;
	const int y_end   = (height + 7) / 8 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width = this->width;
	const int pitch = this->pitch;
	const int hloop = (width + 7) / 8;

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp = luma[1] + y * pitch + 8;
		short* dstp1 = bufy[0] +  y / 2      * pitch + 8;
		short* dstp2 = bufy[1] + (y / 2 + 1) * pitch + 8;

		__asm
		{
			mov			rsi, srcp			// esi = srcp
			mov			rdi, dstp1			// edi = dstp1
			mov			rdx, dstp2			// edx = dstp2
			mov			eax, pitch
			mov			ecx, hloop			// ecx = hloop
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3

align 16
next8columns:
#if !defined(_WIN64)
			push		esi
#else
			mov			r12, rsi
#endif
			movdqa		xmm2, [rsi]
			movdqa		xmm0, [rsi+rax]
			movdqa		xmm1, [rsi+2*rax]
			paddw		xmm2, xmm1
			psraw		xmm2, 1
			psubw		xmm0, xmm2
			add			rsi, rbx

			movdqa		xmm2, [rsi]
			movdqa		xmm3, [rsi+rax]
			movdqa		xmm4, [rsi+2*rax]
			movdqa		xmm5, [rsi+rbx]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdx], xmm2
			movdqa		[rdx+rax], xmm4
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 2
			psraw		xmm2, 2
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi], xmm6
			movdqa		[rdi+rax], xmm7
			lea			rsi, [rsi+4*rax]

			movdqa		xmm0, [rsi]
			movdqa		xmm1, [rsi+rax]
			movdqa		xmm2, [rsi+2*rax]
			movdqa		xmm3, [rsi+rbx]
			movdqa		xmm6, xmm5
			movdqa		xmm7, xmm1
			paddw		xmm5, xmm1
			paddw		xmm1, xmm3
			psraw		xmm5, 1
			psraw		xmm1, 1
			psubw		xmm0, xmm5
			psubw		xmm2, xmm1
			movdqa		[rdx+2*rax], xmm0
			movdqa		[rdx+rbx], xmm2
			paddw		xmm4, xmm0
			paddw		xmm0, xmm2
			psraw		xmm4, 2
			psraw		xmm0, 2
			paddw		xmm6, xmm4
			paddw		xmm7, xmm0
			movdqa		[rdi+2*rax], xmm6
			movdqa		[rdi+rbx], xmm7

#if !defined(_WIN64)
			pop			esi
#else
			mov			rsi, r12
#endif
			add			rsi, 16
			add			rdi, 16
			add			rdx, 16
			sub			ecx, 1
			jnz			next8columns
		}

		// horizontal reflection
		short* p = dstp1;
		for (int i = 0; i < 4; ++i, p += pitch)
			p[-2] = p[2], p[-1] = p[1], p[width] = p[width-2], p[width+1] = p[width-3];
	}

	// vertical reflection
	if (y_start == 0)
		memcpy(bufy[1], bufy[1] + pitch, pitch * sizeof(short));
	if (thread_id == threads - 1 && height % 2 == 0)
		memcpy(bufy[1] + (height / 2 + 1) * pitch, bufy[1] + (height / 2 - 1) * pitch, pitch * sizeof(short));
}

void MosquitoNR::WaveletHorz2(int thread_id)
{
	const int y_start = (height + 15) / 16 *  thread_id      / threads * 8;
	const int y_end   = (height + 15) / 16 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width  = this->width;
	const int pitch  = this->pitch;
	const int hloop1 = (width + 4 + 2 + 3) / 4;
	const int hloop2 = (width + 7) / 8;
	short* work = this->work[thread_id];

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp = bufy[0] + y * pitch + 4;
		short* dstp = bufx[1] + y / 2 * pitch + 8;

		__asm
		{
			// shuffle
			mov			rsi, srcp			// esi = srcp
			mov			rdi, work			// edi = work
			mov			eax, pitch
			mov			ecx, hloop1			// ecx = hloop1
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3
			lea			rdx, [rsi+4*rax]	// edx = srcp + 4 * pitch

align 16
shuffle_next4columns:
			movq		xmm0, qword ptr [rsi]				// 03, 02, 01, 00
			movq		xmm1, qword ptr [rsi+rax]			// 13, 12, 11, 10
			movq		xmm2, qword ptr [rsi+2*rax]			// 23, 22, 21, 20
			movq		xmm3, qword ptr [rsi+rbx]			// 33, 32, 31, 30
			movq		xmm4, qword ptr [rdx]				// 43, 42, 41, 40
			movq		xmm5, qword ptr [rdx+rax]			// 53, 52, 51, 50
			movq		xmm6, qword ptr [rdx+2*rax]			// 63, 62, 61, 60
			movq		xmm7, qword ptr [rdx+rbx]			// 73, 72, 71, 70
			punpcklwd	xmm0, xmm1			// 13, 03, 12, 02, 11, 01, 10, 00
			punpcklwd	xmm2, xmm3			// 33, 23, 32, 22, 31, 21, 30, 20
			punpcklwd	xmm4, xmm5			// 53, 43, 52, 42, 51, 41, 50, 40
			punpcklwd	xmm6, xmm7			// 73, 63, 72, 62, 71, 61, 70, 60
			movdqa		xmm1, xmm0
			movdqa		xmm5, xmm4
			punpckldq	xmm0, xmm2			// 31, 21, 11, 01, 30, 20, 10, 00
			punpckhdq	xmm1, xmm2			// 33, 23, 13, 03, 32, 22, 12, 02
			punpckldq	xmm4, xmm6			// 71, 61, 51, 41, 70, 60, 50, 40
			punpckhdq	xmm5, xmm6			// 73, 63, 53, 43, 72, 62, 52, 42
			movdqa		xmm2, xmm0
			movdqa		xmm3, xmm1
			punpcklqdq	xmm0, xmm4			// 70, 60, 50, 40, 30, 20, 10, 00
			punpckhqdq	xmm2, xmm4			// 71, 61, 51, 41, 31, 21, 11, 01
			punpcklqdq	xmm1, xmm5			// 72, 62, 52, 42, 32, 22, 12, 02
			punpckhqdq	xmm3, xmm5			// 73, 63, 53, 43, 33, 23, 13, 03
			movdqa		[rdi], xmm0
			movdqa		[rdi+16], xmm2
			movdqa		[rdi+32], xmm1
			movdqa		[rdi+48], xmm3
			add			rsi, 8
			add			rdx, 8
			add			rdi, 64
			sub			ecx, 1
			jnz			shuffle_next4columns

			// wavelet transform
			mov			rsi, work
			mov			rdi, dstp			// edi = dstp
			mov			ecx, hloop2			// ecx = hloop2
			add			esi, 64				// esi = work + 32

			movdqa		xmm2, [rsi-32]
			movdqa		xmm0, [rsi-16]
			movdqa		xmm1, [rsi]
			paddw		xmm2, xmm1
			psraw		xmm2, 1
			psubw		xmm0, xmm2
			movdqa		[rdi-16], xmm0

align 16
wavelet_next8columns:
			movdqa		xmm2, [rsi+16]
			movdqa		xmm3, [rsi+32]
			movdqa		xmm4, [rsi+48]
			movdqa		xmm5, [rsi+64]
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdi], xmm2
			movdqa		[rdi+16], xmm4
			movdqa		xmm1, xmm5
			movdqa		xmm2, [rsi+80]
			movdqa		xmm3, [rsi+96]
			movdqa		xmm4, [rsi+112]
			movdqa		xmm5, [rsi+128]
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdi+32], xmm2
			movdqa		[rdi+48], xmm4
			movdqa		xmm1, xmm5
			add			rsi, 128
			add			rdi, 64
			sub			ecx, 1
			jnz			wavelet_next8columns

			// horizontal reflection
			mov			eax, width
			test		eax, 1
			jnz			no_reflection		// if (width % 2 == 0) {
			mov			rdi, dstp			//   edi = dstp
			shl			eax, 3				//   eax = width / 2 * 16
			movdqa		xmm0, [rdi+rax-32]
			movdqa		[rdi+rax], xmm0
no_reflection:								// }
		}
	}
}

void MosquitoNR::WaveletHorz3(int thread_id)
{
	const int y_start = (height + 15) / 16 *  thread_id      / threads * 8;
	const int y_end   = (height + 15) / 16 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width  = this->width;
	const int pitch  = this->pitch;
	const int hloop1 = (width + 4 + 2 + 3) / 4;
	const int hloop2 = (width + 3) / 4;
	short* work = this->work[thread_id];

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp = bufy[0] + y * pitch + 4;
		short* dstp1 = bufx[0] + y / 2 * pitch + 8;
		short* dstp2 = bufx[1] + y / 2 * pitch + 8;

		__asm
		{
			// shuffle
			mov			rsi, srcp			// esi = srcp
			mov			rdi, work			// edi = work
			mov			eax, pitch
			mov			ecx, hloop1			// ecx = hloop1
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3
			lea			rdx, [rsi+4*rax]	// edx = srcp + 4 * pitch

align 16
shuffle_next4columns:
			movq		xmm0, qword ptr [rsi]				// 03, 02, 01, 00
			movq		xmm1, qword ptr [rsi+rax]			// 13, 12, 11, 10
			movq		xmm2, qword ptr [rsi+2*rax]			// 23, 22, 21, 20
			movq		xmm3, qword ptr [rsi+rbx]			// 33, 32, 31, 30
			movq		xmm4, qword ptr [rdx]				// 43, 42, 41, 40
			movq		xmm5, qword ptr [rdx+rax]			// 53, 52, 51, 50
			movq		xmm6, qword ptr [rdx+2*rax]			// 63, 62, 61, 60
			movq		xmm7, qword ptr [rdx+rbx]			// 73, 72, 71, 70
			punpcklwd	xmm0, xmm1			// 13, 03, 12, 02, 11, 01, 10, 00
			punpcklwd	xmm2, xmm3			// 33, 23, 32, 22, 31, 21, 30, 20
			punpcklwd	xmm4, xmm5			// 53, 43, 52, 42, 51, 41, 50, 40
			punpcklwd	xmm6, xmm7			// 73, 63, 72, 62, 71, 61, 70, 60
			movdqa		xmm1, xmm0
			movdqa		xmm5, xmm4
			punpckldq	xmm0, xmm2			// 31, 21, 11, 01, 30, 20, 10, 00
			punpckhdq	xmm1, xmm2			// 33, 23, 13, 03, 32, 22, 12, 02
			punpckldq	xmm4, xmm6			// 71, 61, 51, 41, 70, 60, 50, 40
			punpckhdq	xmm5, xmm6			// 73, 63, 53, 43, 72, 62, 52, 42
			movdqa		xmm2, xmm0
			movdqa		xmm3, xmm1
			punpcklqdq	xmm0, xmm4			// 70, 60, 50, 40, 30, 20, 10, 00
			punpckhqdq	xmm2, xmm4			// 71, 61, 51, 41, 31, 21, 11, 01
			punpcklqdq	xmm1, xmm5			// 72, 62, 52, 42, 32, 22, 12, 02
			punpckhqdq	xmm3, xmm5			// 73, 63, 53, 43, 33, 23, 13, 03
			movdqa		[rdi], xmm0
			movdqa		[rdi+16], xmm2
			movdqa		[rdi+32], xmm1
			movdqa		[rdi+48], xmm3
			add			rsi, 8
			add			rdx, 8
			add			rdi, 64
			sub			ecx, 1
			jnz			shuffle_next4columns

			// wavelet transform
			mov			rsi, work
			mov			rdi, dstp1			// edi = dstp1
			mov			rdx, dstp2			// edx = dstp2
			mov			ecx, hloop2			// ecx = hloop2
			add			rsi, 64				// esi = work + 32

			movdqa		xmm2, [rsi-32]
			movdqa		xmm0, [rsi-16]
			movdqa		xmm1, [rsi]
			paddw		xmm2, xmm1
			psraw		xmm2, 1
			psubw		xmm0, xmm2
			movdqa		[rdx-16], xmm0

align 16
wavelet_next4columns:
			movdqa		xmm2, [rsi+16]
			movdqa		xmm3, [rsi+32]
			movdqa		xmm4, [rsi+48]
			movdqa		xmm5, [rsi+64]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 1
			psraw		xmm3, 1
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdx], xmm2
			movdqa		[rdx+16], xmm4
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 2
			psraw		xmm2, 2
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi], xmm6
			movdqa		[rdi+16], xmm7
			movdqa		xmm0, xmm4
			movdqa		xmm1, xmm5
			add			rsi, 64
			add			rdi, 32
			add			rdx, 32
			sub			ecx, 1
			jnz			wavelet_next4columns

			// horizontal reflection
			mov			eax, width
			test		eax, 1
			jnz			no_reflection		// if (width % 2 == 0) {
			mov			rdi, dstp1			//   edi = dstp1
			mov			rdx, dstp2			//   edx = dstp2
			shl			eax, 3				//   eax = width / 2 * 16
			movdqa		xmm0, [rdi+rax-16]
			movdqa		xmm1, [rdx+rax-32]
			movdqa		[rdi+rax], xmm0
			movdqa		[rdx+rax], xmm1
no_reflection:								// }
		}
	}
}

void MosquitoNR::BlendCoef(int thread_id)
{
	const int y_start = ((height + 15) &~ 15) / 4 *  thread_id      / threads;
	const int y_end   = ((height + 15) &~ 15) / 4 * (thread_id + 1) / threads;
	if (y_start == y_end) return;
	const int pitch = this->pitch;
	const int multiplier = ((128 - restore) << 16) + restore;
	short *dstp = luma[0], *srcp = bufx[0];

	__asm
	{
		mov			rdi, dstp
		mov			rsi, srcp
		mov			eax, pitch
		mov			ebx, y_start
		mov			ecx, y_end
		add			eax, eax			// eax = pitch * sizeof(short)
		sub			ecx, ebx
		imul		ebx, eax			// ebx = y_start * pitch * sizeof(short)
		imul		ecx, eax
		add			rdi, rbx			// edi = luma[0] + y_start * pitch
		add			rsi, rbx			// esi = bufx[0] + y_start * pitch
		shr			ecx, 4				// ecx = (y_end - y_start) * pitch / 8
		movd		xmm6, multiplier
		pshufd		xmm6, xmm6, 0		// xmm6 = [128 - restore, restore] * 4
		mov			edx, 64
		movd		xmm7, edx
		pshufd		xmm7, xmm7, 0		// xmm7 = [64] * 4

align 16
next8pixels:
		movdqa		xmm0, [rdi]			// d7, d6, d5, d4, d3, d2, d1, d0
		movdqa		xmm2, [rsi]			// s7, s6, s5, s4, s3, s2, s1, s0
		movdqa		xmm1, xmm0
		punpcklwd	xmm0, xmm2			// s3, d3, s2, d2, s1, d1, s0, d0
		punpckhwd	xmm1, xmm2			// s7, d7, s6, d6, s5, d5, s4, d4
		pmaddwd		xmm0, xmm6
		pmaddwd		xmm1, xmm6
		paddd		xmm0, xmm7
		paddd		xmm1, xmm7
		psrad		xmm0, 7
		psrad		xmm1, 7
		packssdw	xmm0, xmm1
		movdqa		[rdi], xmm0
		add			rdi, 16
		add			rsi, 16
		sub			ecx, 1
		jnz			next8pixels
	}
}

void MosquitoNR::InvWaveletHorz(int thread_id)
{
	const int y_start = (height + 15) / 16 *  thread_id      / threads * 8;
	const int y_end   = (height + 15) / 16 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int width = this->width;
	const int pitch = this->pitch;
	const int hloop = (width + 3) / 4;
	short* work = this->work[thread_id];

	for (int y = y_start; y < y_end; y += 8)
	{
		short* srcp1 = luma[0] + y / 2 * pitch + 8;
		short* srcp2 = bufx[1] + y / 2 * pitch + 8;
		short* dstp = bufy[0] + y * pitch + 8;

		__asm
		{
			// wavelet transform
			mov			rsi, srcp1			// esi = srcp1
			mov			rdx, srcp2			// edx = srcp2
			mov			rdi, work			// edi = dstp
			mov			ecx, hloop			// ecx = hloop

			movdqa		xmm2, [rdx-16]
			movdqa		xmm0, [rsi]
			movdqa		xmm1, [rdx]
			paddw		xmm2, xmm1
			psraw		xmm2, 2
			psubw		xmm0, xmm2
			movdqa		[rdi], xmm0

align 16
wavelet_next4columns:
			movdqa		xmm2, [rsi+16]
			movdqa		xmm3, [rdx+16]
			movdqa		xmm4, [rsi+32]
			movdqa		xmm5, [rdx+32]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 2
			psraw		xmm3, 2
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdi+32], xmm2
			movdqa		[rdi+64], xmm4
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 1
			psraw		xmm2, 1
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi+16], xmm6
			movdqa		[rdi+48], xmm7
			movdqa		xmm0, xmm4
			movdqa		xmm1, xmm5
			add			rsi, 32
			add			rdx, 32
			add			rdi, 64
			sub			ecx, 1
			jnz			wavelet_next4columns

			// shuffle
			mov			rsi, work			// esi = work
			mov			rdi, dstp			// edi = srcp
			mov			eax, pitch
			mov			ecx, hloop			// ecx = hloop
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3
			lea			rdx, [rdi+4*rax]	// edx = dstp + 4 * pitch

align 16
shuffle_next4columns:
			movdqa		xmm0, [rsi]			// 70, 60, 50, 40, 30, 20, 10, 00
			movdqa		xmm1, [rsi+16]		// 71, 61, 51, 41, 31, 21, 11, 01
			movdqa		xmm2, [rsi+32]		// 72, 62, 52, 42, 32, 22, 12, 02
			movdqa		xmm3, [rsi+48]		// 73, 63, 53, 43, 33, 23, 13, 03
			movdqa		xmm4, xmm0
			movdqa		xmm6, xmm2
			punpcklwd	xmm0, xmm1			// 31, 30, 21, 20, 11, 10, 01, 00
			punpckhwd	xmm4, xmm1			// 71, 70, 61, 60, 51, 50, 41, 40
			punpcklwd	xmm2, xmm3			// 33, 32, 23, 22, 13, 12, 03, 02
			punpckhwd	xmm6, xmm3			// 73, 72, 63, 62, 53, 52, 43, 42
			movdqa		xmm1, xmm0
			movdqa		xmm5, xmm4
			punpckldq	xmm0, xmm2			// 13, 12, 11, 10, 03, 02, 01, 00
			punpckhdq	xmm1, xmm2			// 33, 32, 31, 30, 23, 22, 21, 20
			punpckldq	xmm4, xmm6			// 53, 52, 51, 50, 43, 42, 41, 40
			punpckhdq	xmm5, xmm6			// 73, 72, 71, 70, 63, 62, 61, 60
			movq		qword ptr [rdi], xmm0
			movq		qword ptr [rdi+2*rax], xmm1
			movq		qword ptr [rdx], xmm4
			movq		qword ptr [rdx+2*rax], xmm5
			punpckhqdq	xmm0, xmm0
			punpckhqdq	xmm1, xmm1
			punpckhqdq	xmm4, xmm4
			punpckhqdq	xmm5, xmm5
			movq		qword ptr [rdi+rax], xmm0
			movq		qword ptr [rdi+rbx], xmm1
			movq		qword ptr [rdx+rax], xmm4
			movq		qword ptr [rdx+rbx], xmm5
			add			rsi, 64
			add			rdi, 8
			add			rdx, 8
			sub			ecx, 1
			jnz			shuffle_next4columns
		}
	}

	// vertical reflection
	if (thread_id == threads - 1 && height % 2 == 0)
		memcpy(bufy[0] + height / 2 * pitch, bufy[0] + (height / 2 - 1) * pitch, pitch * sizeof(short));
}

void MosquitoNR::InvWaveletVert(int thread_id)
{
	const int y_start = (height + 7) / 8 *  thread_id      / threads * 8;
	const int y_end   = (height + 7) / 8 * (thread_id + 1) / threads * 8;
	if (y_start == y_end) return;
	const int pitch  = this->pitch;

	for (int y = y_start; y < y_end; y += 8)
	{
		int hloop = (width + 7) / 8;
		short* srcp1 = bufy[0] + y / 2 * pitch + 8;
		short* srcp2 = bufy[1] + y / 2 * pitch + 8;
		short* dstp = luma[1] + (y + 2) * pitch + 8;

		__asm
		{
			mov			rsi, srcp1			// esi = srcp1
			mov			rdx, srcp2			// edx = srcp2
			mov			rdi, dstp			// edi = dstp
			mov			eax, pitch
			add			eax, eax			// eax = pitch * sizeof(short)
			lea			ebx, [eax+2*eax]	// ebx = pitch * sizeof(short) * 3
			lea			ecx, [eax+4*eax]	// ecx = pitch * sizeof(short) * 5

align 16
next8columns:
#if !defined(_WIN64)
			push		edi
#else
			mov			r13, rdi
#endif
			movdqa		xmm2, [rdx]
			movdqa		xmm0, [rsi]
			movdqa		xmm1, [rdx+rax]
			paddw		xmm2, xmm1
			psraw		xmm2, 2
			psubw		xmm0, xmm2
			movdqa		[rdi], xmm0

			movdqa		xmm2, [rsi+rax]
			movdqa		xmm3, [rdx+2*rax]
			movdqa		xmm4, [rsi+2*rax]
			movdqa		xmm5, [rdx+rbx]
			movdqa		xmm6, xmm1
			movdqa		xmm7, xmm3
			paddw		xmm1, xmm3
			paddw		xmm3, xmm5
			psraw		xmm1, 2
			psraw		xmm3, 2
			psubw		xmm2, xmm1
			psubw		xmm4, xmm3
			movdqa		[rdi+2*rax], xmm2
			movdqa		[rdi+4*rax], xmm4
			paddw		xmm0, xmm2
			paddw		xmm2, xmm4
			psraw		xmm0, 1
			psraw		xmm2, 1
			paddw		xmm6, xmm0
			paddw		xmm7, xmm2
			movdqa		[rdi+rax], xmm6
			movdqa		[rdi+rbx], xmm7
			lea			rdi, [rdi+4*rax]

			movdqa		xmm0, [rsi+rbx]
			movdqa		xmm1, [rdx+4*rax]
			movdqa		xmm2, [rsi+4*rax]
			movdqa		xmm3, [rdx+rcx]
			movdqa		xmm6, xmm5
			movdqa		xmm7, xmm1
			paddw		xmm5, xmm1
			paddw		xmm1, xmm3
			psraw		xmm5, 2
			psraw		xmm1, 2
			psubw		xmm0, xmm5
			psubw		xmm2, xmm1
			movdqa		[rdi+2*rax], xmm0
			paddw		xmm4, xmm0
			paddw		xmm0, xmm2
			psraw		xmm4, 1
			psraw		xmm0, 1
			paddw		xmm6, xmm4
			paddw		xmm7, xmm0
			movdqa		[rdi+rax], xmm6
			movdqa		[rdi+rbx], xmm7

#if !defined(_WIN64)
			pop			edi
#else
			mov			rdi, r13
#endif
			add			rsi, 16
			add			rdx, 16
			add			rdi, 16
			sub			hloop, 1
			jnz			next8columns
		}
	}
}
