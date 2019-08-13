/*******************************************************************************
*
*    File Name:  bch_decoder.c
*     Revision:  2.0
*         Date:  March, 2007
*        Email:  nandsupport@micron.com
*      Company:  Micron Technology, Inc.
*
*  Description:  Micron NAND BCH Decoder
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
*                Copyright 2007 Micron Technology, Inc. All rights reserved.
*
* Rev  Author			Date		Changes
* ---  ---------------	----------	-------------------------------
* 1.0  ZS		08/07/2006	Initial release
* 2.0  PF		03/05/2007	Fixed bug that caused some codewords
* 					to not be corrected
* 
* 
/*******************************************************************************/

#include "bch_global.c"

int bb[rr_max] ;	// Syndrome polynomial
int s[rr_max];		// Syndrome values
int syn_error;		// Syndrome error indicator
int count;		// Number of errors
int location[tt_max];	// Error location
int ttx2;		// 2t
int decode_flag;	// Decoding indicator 
	
void parallel_syndrome() {
/* Parallel computation of 2t syndromes.
 * Use the same lookahead matrix T_G_R of parallel computation of parity check bits.
 * The incoming streams are fed into registers from left hand
 */
	int i, j, iii, Temp, bb_temp[rr_max] ;
	int loop_count ;

	// Determine the number of loops required for parallelism.  
	loop_count = ceil(nn_shorten / (double)Parallel) ;
	
	// Serial to parallel data conversion
	for (i = 0; i < Parallel; i++) 
		for (j = 0; j < loop_count; j++) 
			if (i + j * Parallel < nn_shorten)
				data_p[i][j] = recd[i + j * Parallel];
			else
				data_p[i][j] = 0;
	
	// Initialize the parity bits.
	for (i = 0; i < rr; i++)
		bb[i] = 0;
	
	// Compute syndrome polynomial: S(x) = C(x) mod g(x)
	// S(t) = T_G_R S(t-1) + R(t) 
	// Ref: L&C, pp. 225, Fig. 6.11
	for (iii = loop_count - 1; iii >= 0; iii--) {
		for (i = 0; i < rr; i++) {
			Temp = 0;
			for (j = 0; j < rr; j++) 
				if (bb[j] !=0 && T_G_R[i][j] != 0)
					Temp ^= 1 ;
			bb_temp[i] = Temp;
		}
		
		for (i = 0; i < rr; i++)
			bb[i] = bb_temp[i];
		
		for (i = 0; i < Parallel; i++)
			bb[i] = bb[i] ^ data_p[i][iii];
	}
	
	// Computation 2t syndromes based on S(x)
	// Odd syndromes
	syn_error = 0 ;
	for (i = 1; i <= ttx2 - 1; i = i+2) {
	 	s[i] = 0 ;
		for (j = 0; j < rr; j++)
			if (bb[j] != 0)
				s[i] ^= alpha_to[(index_of[bb[j]] + i*j) % nn] ;
		if (s[i] != 0)
			syn_error = 1 ;	// set flag if non-zero syndrome => error
    	}

	// Even syndrome = (Odd syndrome) ** 2
	for (i = 2; i <= ttx2; i = i + 2) {
	 	j = i / 2;
		if (s[j] == 0)
			s[i] = 0;
		else
			s[i] =  alpha_to[(2 * index_of[s[j]]) % nn];
	}
	
	if (Verbose) {
		fprintf(stdout, "# The syndrome from parallel decoder is:\n") ;
		for (i = 1; i <= ttx2; i++)
			fprintf(stdout, "   %4d (%4d) == 0x%04x (0x%x)\n", s[i],index_of[s[i]],s[i], index_of[s[i]]) ;
		fprintf(stdout, "\n\n") ;
	}
}

void decode_bch() {
	register int i, j, elp_sum ;
	int L[ttx2+3];			// Degree of ELP 
	int u_L[ttx2+3];		// Difference between step number and the degree of ELP
	int reg[tt+3];			// Register state
	int elp[ttx2+4][ttx2+4]; 	// Error locator polynomial (ELP)
	int desc[ttx2+4];		// Discrepancy 'mu'th discrepancy
	int u;				// u = 'mu' + 1 and u ranges from -1 to 2*t (see L&C)
	int q;				//

	parallel_syndrome() ;
	
	if (!syn_error) {
		decode_flag = 1 ;	// No errors
		count = 0 ;
	}
	else {	
		// Having errors, begin decoding procedure
		// Simplified Berlekamp-Massey Algorithm for Binary BCH codes
		// 	Ref: Blahut, pp.191, Chapter 7.6 
		// 	Ref: L&C, pp.212, Chapter 6.4
		//
		// Following the terminology of Lin and Costello's book:   
		// 	desc[u] is the 'mu'th discrepancy, where  
		// 	u='mu'+1 and 
		// 	'mu' (the Greek letter!) is the step number ranging 
		// 		from -1 to 2*t (see L&C)
		// 	l[u] is the degree of the elp at that step, and 
		// 	u_L[u] is the difference between the step number 
		// 		and the degree of the elp. 
		
		if (Verbose) fprintf(stdout,"Beginning Berlekamp loop\n");

		// initialise table entries
		for (i = 1; i <= ttx2; i++) 
			s[i] = index_of[s[i]];

		desc[0] = 0;				/* index form */
		desc[1] = s[1];				/* index form */
		elp[0][0] = 1;				/* polynomial form */
		elp[1][0] = 1;				/* polynomial form */
		//elp[2][0] = 1;				/* polynomial form */
		for (i = 1; i < ttx2; i++) {
			elp[0][i] = 0;			/* polynomial form */
			elp[1][i] = 0;			/* polynomial form */
			//elp[2][i] = 0;			/* polynomial form */
		}
		L[0] = 0;
		L[1] = 0;
		//L[2] = 0;
		u_L[0] = -1;
		u_L[1] = 0;
		//u_L[2] = 0;
		u = -1; 
 
		do {
			// even loops always produce no discrepany so they can be skipped
			u = u + 2; 
			if (Verbose) fprintf(stdout,"Loop %d:\n", u);
			if (Verbose) fprintf(stdout,"     desc[%d] = %x\n", u, desc[u]);
			if (desc[u] == -1) {
				L[u + 2] = L[u];
				for (i = 0; i <= L[u]; i++)
					elp[u + 2][i] = elp[u][i]; 
			}
			else {
				// search for words with greatest u_L[q] for which desc[q]!=0 
				q = u - 2;
				if (q<0) q=0;
				// Look for first non-zero desc[q] 
				while ((desc[q] == -1) && (q > 0))
					q=q-2;
				if (q < 0) q = 0;

				// Find q such that desc[u]!=0 and u_L[q] is maximum
				if (q > 0) {
					j = q;
				  	do {
				    		j=j-2;
						if (j < 0) j = 0;
				    		if ((desc[j] != -1) && (u_L[q] < u_L[j]))
				      			q = j;
				  	} while (j > 0);
				}
 
				// store degree of new elp polynomial
				if (L[u] > L[q] + u - q)
					L[u + 2] = L[u];
				else
					L[u + 2] = L[q] + u - q;
 
				// Form new elp(x)
				for (i = 0; i < ttx2; i++) 
					elp[u + 2][i] = 0;
				for (i = 0; i <= L[q]; i++) 
					if (elp[q][i] != 0)
						elp[u + 2][i + u - q] = alpha_to[(desc[u] + nn - desc[q] + index_of[elp[q][i]]) % nn];
				for (i = 0; i <= L[u]; i++) 
					elp[u + 2][i] ^= elp[u][i];

			}
			u_L[u + 2] = u+1 - L[u + 2];
 
			// Form (u+2)th discrepancy.  No discrepancy computed on last iteration 
			if (u < ttx2) {	
				if (s[u + 2] != -1)
					desc[u + 2] = alpha_to[s[u + 2]];
				else 
					desc[u + 2] = 0;

				for (i = 1; i <= L[u + 2]; i++) 
					if ((s[u + 2 - i] != -1) && (elp[u + 2][i] != 0))
			        		desc[u + 2] ^= alpha_to[(s[u + 2 - i] + index_of[elp[u + 2][i]]) % nn];
			 	// put desc[u+2] into index form 
				desc[u + 2] = index_of[desc[u + 2]];	

			}

			if (Verbose) {
				fprintf(stdout,"     deg(elp) = %2d --> elp(%2d):", L[u], u);
				for (i=0; i<=L[u]; i++)
					fprintf(stdout,"  0x%x", elp[u][i]);
				fprintf(stdout,"\n");
				fprintf(stdout,"     deg(elp) = %2d --> elp(%2d):", L[u+2], u+2);
				for (i=0; i<=L[u+2]; i++)
					fprintf(stdout,"  0x%x", elp[u+2][i]);
				fprintf(stdout,"\n");
				fprintf(stdout,"     u_L[%2d] = %2d\n", u, u_L[u]);
				fprintf(stdout,"     u_L[%2d] = %2d\n", u+2, u_L[u+2]);
			}

		} while ((u < (ttx2-1)) && (L[u + 2] <= tt)); 
		if (Verbose) fprintf(stdout,"\n");
		u=u+2;
		L[ttx2-1] = L[u];
		
		if (L[ttx2-1] > tt) 
			decode_flag = 0;
		else {
			// Chien's search to find roots of the error location polynomial
			// Ref: L&C pp.216, Fig.6.1
			if (Verbose) fprintf(stdout,"Chien Search:  L[%d]=%d=%x\n", ttx2-1,L[ttx2-1],L[ttx2-1]);
			if (Verbose) fprintf(stdout,"Sigma(x) = \n");

			if (Verbose) 
				for (i = 0; i <= L[u]; i++) 
					if (elp[u][i] != 0)
						fprintf(stdout,"    %4d (%4d)\n", elp[u][i], index_of[elp[u][i]]);
					else
						fprintf(stdout,"     0\n");

			for (i = 1; i <= L[ttx2-1]; i++) {
				reg[i] = index_of[elp[u][i]];
				if (Verbose) fprintf(stdout,"  reg[%d]=%d=%x\n", i,reg[i],reg[i]);
			}
			count = 0 ;
			// Begin chien search 
			for (i = 1; i <= nn; i++) {
			 	elp_sum = 1 ;
				for (j = 1; j <= L[ttx2-1]; j++) 
					if (reg[j] != -1) {
					 	reg[j] = (reg[j] + j) % nn ;
						elp_sum ^= alpha_to[reg[j]] ;
					}

				// store root and error location number indices
				if (!elp_sum) {
					location[count] = nn - i ;
					if (Verbose) fprintf(stdout,"count: %d location: %d L[ttx2-1] %d\n",
							count,location[count],L[ttx2-1]);
					count++ ;
				}
			}
			
			// Number of roots = degree of elp hence <= tt errors
			if (count == L[ttx2-1]) {   
				decode_flag = 1 ;
				// Correct errors by flipping the error bit
				for (i = 0; i < L[ttx2-1]; i++) 
				 	recd[location[i]] ^= 1 ;	
			}
			// Number of roots != degree of ELP => >tt errors and cannot solve
			else 
				decode_flag = 0 ;
		}
	}
}

int main(int argc,  char** argv)
{	int i, j ;
	int Help ;
	int Input_kk, Output_Syndrome ;			// Input & Output switch
	int in_count, in_v, in_codeword;		// Input statistics
	int decode_success, decode_fail;		// Decoding statistics
	int code_success[kk_max], code_fail[kk_max];	// Decoded and failed words
	int codeword[kk_max], recd_data[kk_max], recd_parity[kk_max] ;
	char in_char;
	
	fprintf(stderr, "# Binary BCH decoder.  Use -h for details.\n\n");
	
	Verbose = 0;
	Input_kk = 0;
	Help = 0;
	mm = df_m;
	tt = df_t;
	Parallel = df_p;
	decode_success = 0; 
	decode_fail = 0;
	for (i=1; i < argc;i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'm': mm = atoi(argv[++i]);
					break;
				case 't': tt = atoi(argv[++i]);
					break;
				case 'p': Parallel = atoi(argv[++i]);
					break;
				case 'k': kk_shorten = atoi(argv[++i]);
					Input_kk = 1;
					if (kk_shorten % 4 != 0) {
						fprintf(stderr, "### k must divide 4.\n\n");
						Help = 1;
					}
					break;
				case 's': Output_Syndrome = 1;
					break;
				case 'v': Verbose = 1;
					break;
				default: Help = 1;
			}
		}
		else 
			Help = 1;
	}
	
	if (Help == 1) {
		fprintf(stdout,"# Usage %s:  BCH decoder\n",argv[0]);
		fprintf(stdout,"    -h:  This help message\n");
		fprintf(stdout,"    -m <field>:  Galois field, GF, for code.  Code length is 2^<field>-1.\n");
		fprintf(stdout,"         The default value is %d for a code length %d.  If the parameter is\n", df_m, (int)pow(2,df_m) - 1);
		fprintf(stdout,"         set to 0, the program will estimate the value based upon the values\n");
		fprintf(stdout,"         chosen for k and t.\n");
		fprintf(stdout,"    -t <correct>:  Correction power of the code.  Default = %d\n", df_t);
		fprintf(stdout,"    -k <data bits>:  Number of data bits to be encoded. Must divide 4.\n");
		fprintf(stdout,"         The default value is the maximum supported by the code which\n");
		fprintf(stdout,"         depends upon the field (-m) and the correction (-t) chosen.\n");
		fprintf(stdout,"    -p <parallel>:  Parallelism in decoder.  Does not effect results but\n");
		fprintf(stdout,"         does change the algorithm used to generate them.  Default = %d\n", df_p);
		fprintf(stdout,"    -s   Syndrome output after the decoded data.  Default disabled. \n");
		fprintf(stdout,"    -v   Verbose mode.  Output detailed information, such as encoded codeword,\n");
		fprintf(stdout,"         received codeword and decoded codeword.  Default disabled. \n");
		fprintf(stdout,"    <stdin>:  character string to decode in hex format.  All other \n");
		fprintf(stdout,"          characters are ignored.  Comments are enclosed in brackets:  { }.\n");
		fprintf(stdout,"          The hex values are converted to binary and taken <data bits> \n");
		fprintf(stdout,"          at a time.\n");
		fprintf(stdout,"    <stdout>:  resulting decoded character string in hex format.\n");
		fprintf(stdout,"    <stderr>:  information about the decode process as well as error messages.\n");
	}
	else {
		nn = (int)pow(2, mm) - 1;
		nn_shorten = nn;
		
		// generate the Galois Field GF(2**mm)
		generate_gf() ;
		
		// Compute the generator polynomial and lookahead matrix for BCH code
		gen_poly() ;
		
		// Check if code is shortened
		if (Input_kk == 1)
			nn_shorten = kk_shorten + rr ;
		else {
			kk_shorten = nn_shorten - rr ;
			// Make the shortened length divide 4
			kk_shorten = kk_shorten - kk_shorten % 4 ;
			nn_shorten = kk_shorten + rr ;
		}
		ttx2 = 2 * tt ;
		
		fprintf(stdout, "{# (m = %d, n = %d, k = %d, t = %d) Binary BCH code.}\n\n", mm, nn_shorten, kk_shorten, tt) ;
		
		// Set input data.	
		in_count = 0;
		in_codeword = 0;
		in_char = getchar();
		while (in_char != EOF) {
			if (in_char=='{') {
				while ((in_char != EOF) && ((char)in_char != '}'))
					in_char = getchar();
			}
			in_v = hextoint(in_char);		
			if (in_v != -1) {
				for (i = 3; i >= 0; i--) {
					if ((int)pow(2,i) & in_v)
						codeword[in_count] = 1 ;
					else
						codeword[in_count] = 0 ;
					in_count++;
				}
			}
			if (in_count == ceil(nn_shorten / (double)4) * 4) {
				in_codeword++ ;
				// Parity check bits
				for (j = kk_shorten; j < nn_shorten; j++)
					recd[j - kk_shorten] = codeword[j] ;
				// Data bits
				for (j = 0; j < kk_shorten; j++)
					recd[j + rr] = codeword[j] ;

				decode_bch() ;
				
				if ( decode_flag == 1 ) {
					decode_success++ ;
					code_success[decode_success] = in_codeword;
					if (count == 0) 
						fprintf(stdout, "{ Codeword %d: No errors.}\n", in_codeword) ;
					else {
						fprintf(stdout, "{ Codeword %d: %d errors found at location:", in_codeword, count) ;
						for (i = count - 1; i >= 0 ; i--)  {
							// Convert error location from systematic form to storage form 
							if (location[i] >= rr)
								location[i] = location[i] - rr;
							else
								location[i] = location[i] + kk_shorten;
							
							fprintf(stdout, " %d", location[i]) ;
						}
						fprintf(stdout, "}");

						printf("\n");
					}
				}
				else {
					decode_fail++ ;
					code_fail[decode_fail] = in_codeword;
					fprintf(stdout, "{ Codeword %d: Unable to decode!}", in_codeword) ;
					printf("\n");
				}
				// Convert decoded data into information data and parity checks
				for (i = 0; i < kk_shorten; i++)
					recd_data[i] = recd[i + rr];
				for (i = 0; i < rr; i++)
					recd_parity[i] = recd[i];
				print_hex_low(kk_shorten, recd_data, stdout);
				if (Output_Syndrome == 1) {
					fprintf(stdout, "    ");
					print_hex_low(rr, recd_parity, stdout);
					if (Verbose) fprintf(stdout,"rr: %d\n",rr);
				}
				fprintf(stdout, "\n\n");
				in_count = 0;
			}
			in_char = getchar();
		}
		
		fprintf(stdout, "{### %d codewords received.}\n", in_codeword) ;
		fprintf(stdout, "{@@@ %d codewords are decoded successfully:}\n{", decode_success) ;
		for (i = 1; i <= decode_success; i++)
			fprintf(stdout, " %d", code_success[i]);
		fprintf(stdout, " }\n");
		fprintf(stdout, "{!!! %d codewords are unable to correct:}\n{", decode_fail) ;
		for (i = 1; i <= decode_fail; i++)
			fprintf(stdout, " %d", code_fail[i]);
		fprintf(stdout, " }\n");
	}
	
	return(0);
}
