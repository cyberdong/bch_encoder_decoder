/*******************************************************************************
*
*    File Name:  error_generator.c
*     Revision:  2.0
*         Date:  March, 2007
*        Email:  nandsupport@micron.com
*      Company:  Micron Technology, Inc.
*
*  Description:  corrupt data streams randomly
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
* 
* Rev  Author		Date		Changes
* ---  ---------------	----------	-------------------------------
* 1.0  ZS		08/07/2006	Initial release
* 2.0  PF		03/05/2007	Added options "-s" and "-r"
* 
* 
******************************************************************************/

#include "bch_global.c"
int codeword[16 * kk_max] ;	// Incoming data

int main(int argc,  char** argv)
{	int i, j ;
	int Error_Number ;	// Number of errors applied
	int Help;
	int rec_mode=0;
	int in_count, in_v;
	int in_count_rec;
	int seed;
	char in_char;
	int wd_cnt[8];
	
	fprintf(stderr, "# Random error generator.  Use -h for details.\n\n");
	
	Help = 0;
	Error_Number = 1;
	seed = 1;
	for (i=1; i < argc;i++) 
	{	if (argv[i][0] == '-') 
		{	switch (argv[i][1]) 
			{	case 'e': Error_Number = atoi(argv[++i]);
					  break;
				case 's': seed = atoi(argv[++i]);
					  break;
				case 'r': rec_mode=1;
			  		break;
				default: Help = 1;
			}
		}
		else 
			Help = 1;
	}
	
	if (Help == 1)
	{	fprintf(stdout,"# Usage %s:  Error generator\n",argv[0]);
		fprintf(stdout,"    -h:  This help message\n");
		fprintf(stdout,"    -r: record based mode, all errors are in a record (terminated by newline).\n");
		fprintf(stdout,"    -e <error>:  Number of errors in codeword.  Default = %d\n", Error_Number);
		fprintf(stdout,"    -s <seed>:  Set the seed for the random number generator.  Default = %d\n", seed);
		fprintf(stdout,"    <stdout>:  resulting corrupted data string in hex format.\n");
		fprintf(stdout,"    <stderr>:  information about the process as well as error messages\n");
	}
	else
	{	in_count = 0;
		in_count_rec = 1;
	  	wd_cnt[0]=wd_cnt[1]=wd_cnt[2]=wd_cnt[3]=wd_cnt[4]=wd_cnt[5]=wd_cnt[6]=wd_cnt[7]=0;
		in_char = getchar();
		srand(seed);
		fprintf(stdout, "{ Seed = %d }\n",seed);
		
		while (in_char != EOF) 
		{	if (in_char=='{') 
			{	while ((in_char != EOF) && ((char)in_char != '}'))
					in_char = getchar();
			}
			in_v = hextoint(in_char);		
			if (in_v != -1)
			{	for (i = 3; i >= 0; i--) 
				{	if ((int)pow(2,i) & in_v)
						codeword[in_count] = 1 ;
					else
						codeword[in_count] = 0 ;
					in_count++;
				}
			} else 
			{	fprintf(stderr,"in_count: %d\n",in_count);
			  	if (rec_mode && in_count)
			    	{ 	fprintf(stdout, "{%6d) %d errors applied.  Errors locations are:}\n{", in_count_rec, Error_Number);
			      		fprintf(stderr, "{%6d) %d errors applied.  Errors locations are:}\n{", in_count_rec, Error_Number);
			      		for (i = 0; i < Error_Number; i++)
					{	j = rand() % in_count ;
				  		codeword[j] = codeword[j] ^ 1;
				  		fprintf(stdout, " %d", j);
				  		fprintf(stderr, " %d", j);
				  		wd_cnt[j*8/in_count]++;
					}
			      		fprintf(stdout, " }\n\n");
			      		fprintf(stderr, " }\n");
			      		print_hex_low(in_count, codeword, stdout);
			      		fprintf(stdout,"\n");
			      		in_count=0;
			      		in_count_rec++;
			    	}
			}
			in_char = getchar();
		}
		if (rec_mode) return(0);
		fprintf(stderr, "# Total number of bits is: %d.\n\n", in_count) ;
		fprintf(stdout, "{%d errors applied.  Error bits locations are:}\n{", Error_Number);
		for (i = 0; i < Error_Number; i++)
		{	j = rand() % in_count ;
			codeword[j] = codeword[j] ^ 1;
			fprintf(stdout, " %d", j);
			wd_cnt[j*8/in_count]++;
		}
		fprintf(stdout, " }\n\n");
		fprintf(stderr, "word counts:");
		for (i=0;i<8;i++)
			fprintf(stderr," %d",wd_cnt[i]);
		fprintf(stderr,"\n");
		print_hex_low(in_count, codeword, stdout);
	}
	
	return(0);
}
