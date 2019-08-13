/*******************************************************************************
*
*    File Name:  bch_global.c
*     Revision:  1.0
*         Date:  August, 2006
*        Email:  nandsupport@micron.com
*      Company:  Micron Technology, Inc.
*
*  Description:  Micron NAND BCH Global Package
*
*     Function:   1. Create Galois Field
*		  2. Create Generator Polynomial
*		  3. Create Parallel Generator Polynomial 
*
*   References: 
* 		  1. Error Control Coding, Lin & Costello, 2nd Ed., 2004
* 		  2. Error Control Codes, Blahut, 1983
* 		  3. Parallel CRC, Shieh, 2001
*
**
*   Disclaimer   This software code and all associated documentation, comments or other 
*  of Warranty:  information (collectively "Software") is provided "AS IS" without 
*                warranty of any kind. MICRON TECHNOLOGY, INC. ("MTI") EXPRESSLY 
*                DISCLAIMS ALL WARRANTIES EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
*                TO, NONINFRINGEMENT OF THIRD PARTY RIGHTS, AND ANY IMPLIED WARRANTIES 
*                OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE. MTI DOES NOT 
*                WARRANT THAT THE SOFTWARE WILL MEET YOUR REQUIREMENTS, OR THAT THE 
*                OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE. 
*                FURTHERMORE, MTI DOES NOT MAKE ANY REPRESENTATIONS REGARDING THE USE OR 
*                THE RESULTS OF THE USE OF THE SOFTWARE IN TERMS OF ITS CORRECTNESS, 
*                ACCURACY, RELIABILITY, OR OTHERWISE. THE ENTIRE RISK ARISING OUT OF USE 
*                OR PERFORMANCE OF THE SOFTWARE REMAINS WITH YOU. IN NO EVENT SHALL MTI, 
*                ITS AFFILIATED COMPANIES OR THEIR SUPPLIERS BE LIABLE FOR ANY DIRECT, 
*                INDIRECT, CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES (INCLUDING, 
*                WITHOUT LIMITATION, DAMAGES FOR LOSS OF PROFITS, BUSINESS INTERRUPTION, 
*                OR LOSS OF INFORMATION) ARISING OUT OF YOUR USE OF OR INABILITY TO USE 
*                THE SOFTWARE, EVEN IF MTI HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH 
*                DAMAGES. Because some jurisdictions prohibit the exclusion or 
*                limitation of liability for consequential or incidental damages, the 
*                above limitation may not apply to you.
*
*                Copyright 2006 Micron Technology, Inc. All rights reserved.
*
*
* Rev  Author		Date		Changes
* ---  ---------------	----------	-------------------------------
* 1.0  ZS		08/07/2006	Initial release
* 2.0  PF		03/05/2007	Expanded constants to allow 
* 					larger fields
* 
* 
/*******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define mm_max  15         	/* Dimension of Galoise Field */
#define nn_max  32768        	/* Length of codeword, n = 2**m - 1 */
#define tt_max  20          	/* Number of errors that can be corrected */
#define kk_max  32768        	/* Length of information bit, kk = nn - rr  */
#define rr_max  1000		/* Number of parity checks, rr = deg[g(x)] */
#define parallel_max  32	/* Number of parallel encoding/syndrome computations */
#define DEBUG  0

/* Default values */
int df_m = 13;              	// BCH code over GF(2**mm)
int df_t = 4;              	// Number of errors that can be corrected
int df_p = 8;              	// Number of substreams to calculate in parallel

int mm, nn, kk, tt, rr;		// BCH code parameters
int nn_shorten, kk_shorten;	// Shortened BCH code
int Parallel ;			// Parallel processing
int Verbose ;			// Mode indicator
int p[mm_max + 1], alpha_to[nn_max], index_of[nn_max] ;	// Galois field
int gg[rr_max] ;		// Generator polynomial
int T_G[rr_max][rr_max], T_G_R[rr_max][rr_max];		// Parallel lookahead table
int T_G_R_Temp[rr_max][rr_max] ; 
int data[kk_max], data_p[parallel_max][kk_max], recd[nn_max] ;	// Information data and received data

int hextoint(char hex)
// Convert HEX number to Integer
{	int r, h;
	r = -1;
	h = (int)hex;
	if ((h >= 97) && (h <= 102)) 
		r = h - 87;
	else if ((h >= 65) && (h <= 70)) 
		r = h - 55;
	else if ((h >= 48) && (h <= 57)) 
		r = h - 48;
	return r;
}

char inttohex(int i)
// Convert Integer number to HEX
{	char r;
	if (i > 9)
		r = (char)(55 + i);
	else
		r = (char)(48 + i);
	return r;
}

void print_hex(int length, int Binary_data[length], FILE *std)
// Print the binary data in HEX form
// 1100 1010 = 5 3
{	int i, j, l, v;
	l = ceil((double)length / 4);
	
	for (j = l - 1; j >= 0; j--) 
	{	v = 0;
		for(i = 3; i >= 0; i--)
			v = v + (int)Binary_data[j * 4 + i] * (int)pow(2, i);
		fprintf(std, " %c", inttohex(v));
	}
}

void print_hex_low(int length, int Binary_data[length], FILE *std)
// Print the binary data in HEX form from low to high order
// 1100 1010 = C A
{	int i, j, l, v;
	l = ceil((double)length / 4);
	
	for (j = 0; j < l; j++) 
	{	v = 0;
		for(i = 0; i <= 3; i++)
			v = v + (int)Binary_data[j * 4 + i] * (int)pow(2, 3-i);
		fprintf(std, "%c", inttohex(v));
	}
}

void generate_gf()
/* Generate GF(2**mm) from the primitive polynomial p(X) in p[0]..p[mm]
   The lookup table looks like:  
   index -> polynomial form:   alpha_to[ ] contains j = alpha**i;
   polynomial form -> index form:  index_of[j = alpha**i] = i
   alpha_to[1] = 2 is the primitive element of GF(2**mm)
 */
{	int i;
	int mask ;	// Register states
	
	// Primitive polynomials
   	for (i = 1; i < mm; i++)
		p[i] = 0;
	p[0] = p[mm] = 1;
	if (mm == 2)		p[1] = 1;
	else if (mm == 3)	p[1] = 1;
	else if (mm == 4)	p[1] = 1;
	else if (mm == 5)	p[2] = 1;
	else if (mm == 6)	p[1] = 1;
	else if (mm == 7)	p[1] = 1;
	else if (mm == 8)	p[4] = p[5] = p[6] = 1;
	else if (mm == 9)	p[4] = 1;
	else if (mm == 10)	p[3] = 1;
	else if (mm == 11)	p[2] = 1;
	else if (mm == 12)	p[3] = p[4] = p[7] = 1;
	else if (mm == 13)	p[1] = p[2] = p[3] = p[5] = p[7] = p[8] = p[10] = 1;	// 25AF
	// else if (mm == 13)	p[1] = p[3] = p[4] = 1;
	else if (mm == 14)	p[2] = p[4] = p[6] = p[7] = p[8] = 1;	// 41D5
	// else if (mm == 14)	p[1] = p[11] = p[12] = 1;
	else if (mm == 15)	p[1] = 1;
	else if (mm == 16)	p[2] = p[3] = p[5] = 1;
	else if (mm == 17)	p[3] = 1;
	else if (mm == 18)	p[7] = 1;
	else if (mm == 19)	p[1] = p[5] = p[6] = 1;
	else if (mm == 20)	p[3] = 1;
	
	if (Verbose)
	{	fprintf(stderr, "# The Galois field is GF(2**%d);\n\n", mm);
		fprintf(stderr, "# The primitive polynomial is: p(x) = ");
		for (i = 0; i <= mm; i++) 
		{	fprintf(stderr, " %d", p[i]);
		}
		fprintf(stderr, "\n\n");
	}
	
	// Galois field implementation with shift registers
	// Ref: L&C, Chapter 6.7, pp. 217
	mask = 1 ;
	alpha_to[mm] = 0 ;
	for (i = 0; i < mm; i++)
	{ 	alpha_to[i] = mask ;
		index_of[alpha_to[i]] = i ;
		if (p[i] != 0)
			alpha_to[mm] ^= mask ;
		
		mask <<= 1 ;
	}
	
	index_of[alpha_to[mm]] = mm ;
	mask >>= 1 ;
	for (i = mm + 1; i < nn; i++)
	{ 	if (alpha_to[i-1] >= mask)
			alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i-1] ^ mask) << 1) ;
		else alpha_to[i] = alpha_to[i-1] << 1 ;
		
		index_of[alpha_to[i]] = i ;
	}
	index_of[0] = -1 ;
	
	// Print out the Galois Field
	if (Verbose)
	{	fprintf(stderr, "# Look-up tables for GF(2**%2d)\n", mm) ;
		fprintf(stderr, "  i   alpha_to[i]  index_of[i]\n") ;
		for (i=0; i<=nn; i++)
			fprintf(stderr, "%3d      %3d          %3d\n", i, alpha_to[i], index_of[i]) ;
		fprintf(stderr, "\n") ;
	}
}


void gen_poly()
/* Compute generator polynomial of the tt-error correcting Binary BCH code 
 * g(x) = LCM{M_1(x), M_2(x), ..., M_2t(x)},
 * where M_i(x) is the minimal polynomial of alpha^i by cyclotomic cosets
 */
{	int gen_roots[nn + 1], gen_roots_true[nn + 1] ; 	// Roots of generator polynomial
	int i, j, iii, jjj, Temp ;
		
	// Initialization of gen_roots
	for (i = 0; i <= nn; i++) 
	{	gen_roots_true[i] = 0;
		gen_roots[i] = 0;
	}

	// Cyclotomic cosets of gen_roots
   	for (i = 1; i <= 2*tt ; i++)
	{	for (j = 0; j < mm; j++) 
		{	Temp = ((int)pow(2, j) * i) % nn;
			gen_roots_true[Temp] = 1;
		}
	}
	
   	rr = 0;		// Count the number of parity check bits
   	for (i = 0; i < nn; i++) 
	{	if (gen_roots_true[i] == 1) 
		{	rr++;
			gen_roots[rr] = i;
		}
	}
	kk = nn - rr;
	
	// Compute generator polynomial based on its roots
	gg[0] = 2 ;	// g(x) = (X + alpha) initially
	gg[1] = 1 ;
	for (i = 2; i <= rr; i++)
	{ 	gg[i] = 1 ;
		for (j = i - 1; j > 0; j--)
		if (gg[j] != 0)  
			gg[j] = gg[j-1]^ alpha_to[(index_of[gg[j]] + index_of[alpha_to[gen_roots[i]]]) % nn] ;
		else 
			gg[j] = gg[j-1] ;
		gg[0] = alpha_to[(index_of[gg[0]] + index_of[alpha_to[gen_roots[i]]]) % nn] ;
	}
	
	if (Verbose)
	{	fprintf(stderr, "# The Generator Polynomial is:\n") ;
		for (i=0; i <= rr; i++)  
			fprintf(stderr, " %d", gg[i]) ;
		fprintf(stderr, "\n\n") ;
	}
	
	// for parallel encoding and syndrome computation
	// Max parallalism is rr
	if (Parallel > rr)
		Parallel = rr ;
	
	// Construct parallel lookahead matrix T_g, and T_g**r from gg(x)
	// Ref: Parallel CRC, Shieh, 2001
	for (i = 0; i < rr; i++)
	{	for (j = 0; j < rr; j++)
			T_G[i][j] = 0;
	}
	
	for (i = 1; i < rr; i++)
		T_G[i][i-1] = 1 ;

	for (i = 0; i < rr; i++)
		T_G[i][rr - 1] = gg[i] ;
	
	for (i = 0; i < rr; i++)
	{	for (j = 0; j < rr; j++)
			T_G_R[i][j] = T_G[i][j];
	}
	
	// Compute T_G**R Matrix
	for (iii = 1; iii < Parallel; iii++)
	{	for (i = 0; i < rr; i++)
		{	for (j = 0; j < rr; j++)
			{	Temp = 0;
				for (jjj = 0; jjj < rr; jjj++)
					Temp = Temp ^ T_G_R[i][jjj] * T_G[jjj][j];
				
				T_G_R_Temp[i][j] = Temp;
			}
		}
		
		for (i = 0; i < rr; i++)
		{	for (j = 0; j < rr; j++)
				T_G_R[i][j] = T_G_R_Temp[i][j];
		}
	}
}
