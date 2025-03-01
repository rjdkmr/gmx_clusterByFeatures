/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014-2025  Rajendra Kumar
 *
 * gmx_clusterByFeatures uses hole program for which documentation is available in the
 * following link:
 * http://www.csb.yale.edu/userguides/graphics/hole/doc/hole_d00.html
 * Please cite the original publication of hole:
 * O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)
 * The Pore Dimensions of Gramicidin A
 * Biophysical Journal 65:2455-2460.
 *
 * gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gmx_clusterByFeatures is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include "parseData.h"

void freeCharArray(char **a, int m) {
    int i;
    for (i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
}

int* extract_coulmn_integer(char *str, int col_min, int col_max)	{
	/*
	 Similar to extract_coulmn_double but extract integer data.
	 */

	int *data=NULL;
	char *buffer=NULL, **str_data=NULL;
	int size = (col_max-col_min)+1;
	int i=0,n=0, num;

	buffer = strdup(str);
	remove_leading_white_space(buffer);
	str_data = split_by_space(buffer, &num);

	data = (int *) malloc (sizeof(int)*size);
	for(i=(col_min-1);i<col_max;i++){
		data[n] = atoi(str_data[i]);
		n++;
	}
	
	freeCharArray(str_data, num);
	return data;
}




double* extract_coulmn_double(char *str, int col_min, int col_max)	{
	/*
	 return pointer to array of size [(col_max-col_min)+1]
	 array contains floating point numbers range from from col_min to col_max of the given line/string.
	 For example;
	 ---------------------------------------
	 	double *out;
	 	char *line = "    24 U-A       0.07      0.02      0.49     -2.48     -7.63     -0.38";
	 	char *input=NULL;
	 	input = strdup(line);
		out = extract_coulmn_double(input, 3, 5);
		printf("%15.3f%15.3f%15.3f\n",out[0], out[1], out[2]);
	-----------------------------------------
		OUTPUT $          0.070          0.020          0.490
	 */

	double *data=NULL;
	char *buffer=NULL, **str_data=NULL;
	int size = (col_max-col_min)+1;
	int i=0,n=0, num;

	buffer = strdup(str);
	remove_leading_white_space(buffer);
	str_data = split_by_space(buffer, &num);

	data = (double *) malloc (sizeof(double)*size);
	for(i=(col_min-1);i<col_max;i++){
		data[n] = strtof(str_data[i],NULL);
		n++;
	}
	free(buffer);
	free(str_data);
	return data;
}

cbool is_first_numerics(char *str)	{
	char *buffer;
	cbool inumber=FALSE;
	buffer = strdup(str);
	remove_leading_white_space(buffer);
	if(isdigit(buffer[0]))
		inumber=TRUE;
	if((buffer[0]=='-') && (isdigit(buffer[1])))
			inumber=TRUE;
	free(buffer);
	return inumber;
}

char* extract_digits(char *str)	{
	/*
	 Extract integer embedded in the given word string.
	 After extracting the digits, one can directly convert it to integer type by atoi() function
	 For example:
	 -------------------------------
	 char **split_data=NULL;
	 char *line = "10   (0.034) ...1>-:..51_:[.RC]C-----G[.RG]:..63_:-<...1 (0.040)     |"
	 char *input = NULL:
	 char *only_digit=NULL;

	 input = strdup(line);
	 split_data = split_by_char(input, ":");				//Splitting with semi-colon.
	 only_digits = extract_digits(split_data[1]);			//Extracting digits from "..51_"
	 printf("String Type: %s\n",only_digits);
	 printf("Integer Type: %d\n",atoi(only_digits));
	 ---------------------------------
	 OUTPUT $   51

	 */

	char *final=NULL;
	int i = 0, n=0;
	final = (char *) malloc (sizeof(char));
	while(1)	{

		if(str[i]=='\0')	{
			final  = (char*) realloc (final, (sizeof(char)*(n+1)));
			final[n] = '\0';
			break;
		}

		if((i==0) && isdigit(str[i]))	{
			final[n] = str[i];
			n++;
		}

		if((i>0) && isdigit(str[i]))	{
			final  = (char*) realloc (final, (sizeof(char)*(n+1)));
			final[n] = str[i];
			n++;
		}

		i++;
	}
	return final;
}

char** split_by_char(char *OrigStr, char *delimeter)	{
	/*
	 Similar to split sub-routine in perl
	 Split input line with the given character
	 and return the list of words in the given line

	 For example:
	 -------------------------------
	 char **split_data=NULL;
	 char *line = "10   (0.034) ...1>-:..51_:[.RC]C-----G[.RG]:..63_:-<...1 (0.040)     |";
	 char *input=NULL;
	 input = strdup(line);
	 char *only_digits=NULL;

	 split_data = split_by_char(input, ":");				//Splitting with semi-colon.
	 printf("%s ## %s\n",split_data[1], split_data[3]);
	 only_digits = extract_digits(split_data[1]);			//Extracting digits from "..51_"
	 printf("String Type: %s\n",only_digits);
	 printf("Integer Type: %d\n",atoi(only_digits));
	 ---------------------------------
	 OUTPUT $..51_ ## ..63_
	 OUTPUT $String Type: 51
	 OUTPUT $Integer Type: 51
	 */

	char *str = NULL;
	char **final = NULL;
	char *buffer = NULL;
	int n = 0;

	if(OrigStr!=NULL)	{
		str = (char *) malloc(sizeof(char) * strlen(OrigStr)+1);
		strcpy(str,OrigStr);

		final = (char **) malloc (sizeof(char*));
		buffer = strtok(str,delimeter);
		final[n] = strdup(buffer);
		n++;

		while(1)	{
			buffer = strtok(NULL, delimeter);

			if(buffer==NULL)
				break;

			final  = (char**) realloc (final, (sizeof(char*)*(n+1)));
			final[n] = strdup(buffer);
			n++;
		}
	}
	
	free(str);
    return final;
}


char** split_by_space(char *OrigStr, int *num)	{
	/*
	 Split input line with the white spaces
	 and return the list of words in the given line
	 */

	char *str = NULL;
	char **final = NULL;
	char *buffer = NULL;
	int n = 0;

	if(OrigStr!=NULL)	{
		str = (char *) malloc(sizeof(char) * strlen(OrigStr)+1);
		strcpy(str,OrigStr);

		final = (char **) malloc (sizeof(char*));
		buffer = strtok(str," \t\n\v\f\r");
		final[n] = strdup(buffer);
		n++;

		while(1)	{
			buffer = strtok(NULL," \t\n\v\f\r");

			if(buffer==NULL)
				break;

			final  = (char**) realloc (final, (sizeof(char*)*(n+1)));
			final[n] = strdup(buffer);
			n++;
		}
	}
	*num = n;
    free(str);
	return final;
}

void remove_leading_white_space(char *str)	{
	/*
	 Remove any leading white space from the input line/string
	 and replace input string with new line/string.
	 */

	char *final=NULL;
	int i = 0, n =0;
	cbool got_char = FALSE;


	final = (char *) malloc (sizeof(char));
	while(1)	{

		if (str[i]=='\0')	{
			final  = (char*) realloc (final, (sizeof(char)*(n+1)));
			final[n] = '\0';
			break;
		}

		if ((!got_char) && (!isspace(str[i])))	{
				got_char = TRUE;
		}

		if(got_char) {
			if (n==0)
				final[n] = str[i];

			if(n!=0)	{
				final  = (char*) realloc (final, (sizeof(char)*(n+1)));
				final[n] = str[i];
			}
			n++;
		}
		i++;
	}
	strcpy(str,final);
	free(final);
}

char* get_line(FILE *fp)	{
	/*
	 Return one line from the file
	 This function is used in get_all_lines.
	 Each line/string end with '\0' null terminator.
	 */

	char c, *str =NULL;
	char *final=NULL;
	int count = 0;

	str = (char *) malloc (sizeof(char));

	do {

		//Getting character from file
		c = fgetc(fp);

		if(c==EOF)	{
			str  = (char*) realloc (str, (sizeof(char)*(count+1)));
			str[count] = '\0';
			break;
		}

		//Checking end of the line
		if((c=='\n') && (count!=0))	{
			str  = (char*) realloc (str, (sizeof(char)*(count+1)));
			str[count] = c;
			str  = (char*) realloc (str, (sizeof(char)*(count+2)));
			str[count+1] = '\0';
			break;
		}

		//First character, no memory allocation
		if (count==0)	{
				str[count] = c;

				//If first character either new line or end of file
				if((c=='\n') || feof(fp))	{
					str  = (char*) realloc (str, (sizeof(char)*(count+2)));
					str[count+1] = '\0';
				}

		}
		//Any other characters
		else	{
			str  = (char*) realloc (str, (sizeof(char)*(count+1)));
			str[count] = c;
		}

		count++;
	} while (c != EOF);

	final= strdup(str);
	free(str);
	return final;
}

char** get_all_lines(FILE *fp, int *num_line)	{
	/*
	 Return lines from the input file stream.
	 Each line/string end with '\0' null terminator.
	 num_line is output and store the total number of line in the given file.
	 For example:
	 -----------------------------------------
	 	FILE *fp=NULL;
		char **lines=NULL;
		int num_line = 0, i =0;
		fp = fopen("any_file.dat", "r");
		lines = get_all_lines(fp, &num_line);			//Pointer of num_line is passed

		printf("Total number of lines in file: %d",num_line);
		for (i=0;i<num_line;i++){
			printf("%s",lines[i]);
		}
	-------------------------------------------
	 */

	char *buffer=NULL;
	int count = 0;
	char **data = 0;

	data = (char **) malloc (sizeof(char*));
	do {
			buffer = get_line(fp);

			if (buffer==NULL)
				continue;

			if(count==0)	{
				data[count] = buffer;

			}
			else	{
				data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
				data[count] = buffer;
			}

			count++;
	} while(!feof(fp));

	*num_line = (count-1);
	return data;
}

char** get_block_lines(FILE *fp, char *delimeter, int *num_line)	{

	char *buffer=NULL;
	int count = 0;
	char **data = 0;

	data = (char **) malloc (sizeof(char*));
	while(1) {
			buffer = get_line(fp);

			if (buffer==NULL)
				continue;

			if(strstr(buffer,delimeter)!=NULL)
				break;

			if(count==0)	{
				data[count] = buffer;

			}
			else	{
				data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
				data[count] = buffer;
			}

			count++;
	}

	if(buffer!=NULL){
		data  = (char**) realloc (data, (sizeof(char*)*(count+1)));
		data[count] = buffer;
		count++;
	}

	*num_line = (count-1);
	return data;
}
