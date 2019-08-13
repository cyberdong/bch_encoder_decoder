/*******************************************************************************
*
*    File Name:  bch_encoder.c
*     Revision:  1.0
*         Date:  August, 2006
*        Email:  nandsupport@micron.com
*      Company:  Micron Technology, Inc.
*
*  Description:  Micron NAND BCH Encoder
*
*   References: 
* 		  1. Error Control Coding, Lin & Costello, 2nd Ed., 2004
* 		  2. Error Control Codes, Blahut, 1983
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
* Rev  Author			Date		Changes
* ---  ---------------	----------	-------------------------------
* 1.0  ZS		08/07/2006	Initial release
* 
* 
/*******************************************************************************/


#include "bch_global.c"

int bb[rr_max] ;		// Parity checks

void parallel_encode_bch()
/* Parallel computation of n - k parity check bits.
 * Use lookahead matrix T_G_R.
 * The incoming streams are fed into registers from the right hand
 */
{	int i, j, iii, Temp, bb_temp[rr_max] ;
	int loop_count ;
	
	// Determine the number of loops required for parallelism.  
	loop_count = ceil(kk_shorten / (double)Parallel) ;	
	
	// Serial to parallel data conversion
	for (i = 0; i < Parallel; i++)
	{	for (j = 0; j < loop_count; j++)
		{	if (i + j * Parallel < kk_shorten)
				data_p[i][j] = data[i + j * Parallel];
			else
				data_p[i][j] = 0;
		}
	}
	
	// Initialize the parity bits.
	for (i = 0; i < rr; i++)
		bb[i] = 0;
	
	// Compute parity checks
	// S(t) = T_G_R [ S(t-1) + M(t) ]
	// Ref: Parallel CRC, Shieh, 2001
	for (iii = loop_count - 1; iii >= 0; iii--)
	{	for (i = 0; i < rr; i++)
			bb_temp[i] = bb[i] ;
		for (i = Parallel - 1; i >= 0; i--)
			bb_temp[rr - Parallel + i] = bb_temp[rr - Parallel + i] ^ data_p[i][iii];
		
		for (i = 0; i < rr; i++)
		{	Temp = 0;
			for (j = 0; j < rr; j++)
				Temp = Temp ^ (bb_temp[j] * T_G_R[i][j]);
			bb[i] = Temp;
		}
	}
	
}

int main(int argc,  char** argv)
{	int i ;
	int Help ;
	int Input_kk ;				// Input indicator
	int in_count, in_v, in_codeword;	// Input statistics
	char in_char;
	
	fprintf(stderr, "# Binary BCH encoder.  Use -h for details.\n\n");
	
	Verbose = 0;
	Input_kk = 0;
	Help = 0;
	mm = df_m;
	tt = df_t;
	Parallel = df_p;
	for (i = 1; i < argc;i++) 
	{	if (argv[i][0] == '-') 
		{	switch (argv[i][1]) 
			{	case 'm': mm = atoi(argv[++i]);
					if (mm > mm_max)
						Help = 1;
					break;
				case 't': tt = atoi(argv[++i]);
					break;
				case 'p': Parallel = atoi(argv[++i]);
					break;
				case 'k': kk_shorten = atoi(argv[++i]);
					if (kk_shorten % 4 != 0)
					{	fprintf(stderr, "### k must divide 4.\n\n");
						Help = 1;
					}
					Input_kk = 1;
					break;
				case 'v': Verbose = 1;
					break;
				default: Help = 1;
			}
		}
		else 
			Help = 1;
	}
	
	if (Help == 1)
	{	fprintf(stdout,"# Usage %s:  BCH encoder\n", argv[0]);
		fprintf(stdout,"    -h:  This help message\n");
		fprintf(stdout,"    -m <field>:  Galois field, GF, for code.  Code length is 2^<field>-1.\n");
		fprintf(stdout,"         The default value is %d for a code length %d.  If the parameter is\n", df_m, (int)pow(2,df_m) - 1);
		fprintf(stdout,"         set to 0, the program will estimate the value based upon the values\n");
		fprintf(stdout,"         chosen for k and t.\n");
		fprintf(stdout,"    -t <correct>:  Correction power of the code.  Default = %d\n",df_t);
		fprintf(stdout,"    -k <data bits>:  Number of data bits to be encoded. Must divide 4.\n");
		fprintf(stdout,"         The default value is the maximum supported by the code which\n");
		fprintf(stdout,"         depends upon the field (-m) and the correction (-t) chosen.\n");
		fprintf(stdout,"    -p <parallel>:  Parallelism in encoder.  Does not effect results but\n");
		fprintf(stdout,"         does change the algorithm used to generate them.  Default = %d\n", df_p);
		fprintf(stdout,"    -v   Verbose mode.  Output detailed information, such as encoded codeword,\n");
		fprintf(stdout,"         received codeword and decoded codeword.  Default disabled. \n");
		fprintf(stdout,"    <stdin>:  character string to encode in hex format.  All other \n");
		fprintf(stdout,"          characters are ignored.  Comments are enclosed in brackets:  { }.\n");
		fprintf(stdout,"          The hex values are converted to binary and taken <data bits> \n");
		fprintf(stdout,"          at a time.\n");
		fprintf(stdout,"    <stdout>:  resulting encoded character string in hex format.\n");
		fprintf(stdout,"    <stderr>:  information about the encode process as well as error messages.\n");
	}
	else
	{	nn = (int)pow(2, mm) - 1 ;
		nn_shorten = nn ;
		
		// generate the Galois Field GF(2**mm)
		generate_gf() ;
		
		// Compute the generator polynomial and lookahead matrix for BCH code
		gen_poly() ;
		
		// Check if code is shortened
		if (Input_kk == 1)
			nn_shorten = kk_shorten + rr ;
		else
		{	kk_shorten = nn_shorten - rr ;
			// Make the shortened length divide 4
			kk_shorten = kk_shorten - kk_shorten % 4 ;
			nn_shorten = kk_shorten + rr ;
		}
		
		fprintf(stdout, "{# (m = %d, n = %d, k = %d, t = %d, r = %d) Binary BCH code.}\n", mm, nn_shorten, kk_shorten, tt, rr) ;
		
		// Read in data stream
		in_count = 0;
		in_codeword = 0;
		
		in_char = getchar();
		while (in_char != EOF) 
		{	if (in_char=='{') 
			{	while ((in_char != EOF) && ((char)in_char != '}'))
					in_char = getchar();
			}
			in_v = hextoint(in_char);		
			if (in_v != -1)
			{	for (i = 3; i >= 0; i--) 
				{	if ((int)pow(2,i) & in_v)
						data[in_count] = 1 ;
					else
						data[in_count] = 0 ;
					
					in_count++;
				}
			}
			if (in_count == kk_shorten) 
			{	in_codeword++ ;
				
				parallel_encode_bch() ;
				
				print_hex_low(kk_shorten, data, stdout);
				fprintf(stdout, "    ");
				print_hex_low(rr, bb, stdout);
				fprintf(stdout, "\n") ;
				
				in_count = 0;
			}
			in_char = getchar();
		
			// For last codeword
			if (in_char == EOF && in_count > 0) 
			{	in_codeword++ ;
				// Pad zeros
				for (i = in_count; i < kk_shorten; i++)
					data[i] = 0;
				
				parallel_encode_bch() ;
				
				print_hex_low(kk_shorten, data, stdout);
				fprintf(stdout, "    ");
				print_hex_low(rr, bb, stdout);
				fprintf(stdout, "\n") ;
				in_count = 0;
			}
		}
		fprintf(stdout, "\n{### %d words encoded.}\n", in_codeword) ;
	}
	
	return(0);
}
