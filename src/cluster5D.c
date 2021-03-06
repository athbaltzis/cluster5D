#include <stdio.h>
#include <math.h>
#include <dirent.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#define		YES		1

#define		NO		0

#define 	DIMENSIONS		40

#define		MIN_ADD		1000000

#define		MAX_SD_POINTER		40


/****************************************/
/*					*/
/*Variable declarations			*/
/*					*/
/****************************************/


clock_t	begin, end;
double	time_spent;

int AUTO_FILENAME = NO;
int	AUTO_FRACTION = YES;
int	FRACTION = 0;
int VERBOSE = NO;

int original_matrix[DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1];
int	matrix[DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1][DIMENSIONS+1];
int projection[DIMENSIONS+1][DIMENSIONS+1];
int DIM;
float cutoff;
void recursion (int, int, int, int, int);
void smoothing (int, int, int, int, int);
int add = MIN_ADD;
int main(int argc,char *argv[])
{

int	count , count_dpca , count_pca ;
DIR	*dir;
struct dirent	*dp;
char	line[10000];

int	num;
FILE	*fp;
float	pca1,pca2,pca3,pca4,pca5;
float	pca1_max, pca1_min;
float	pca2_max, pca2_min;
float	pca3_max, pca3_min;
float	pca4_max, pca4_min;
float	pca5_max, pca5_min;
int	have_limits;

int	pointer1;
int	pointer2;
int	pointer3;
int	pointer4;
int	pointer5;

int	i, k, j, l, m;
int value,value1;

float	sum;
float	N;
float	mean;
float	val;
float	sum_sd;
float	variance , sd;

int	frames_count;	

double	sd_pointer;
float	local_val;
float	local_sum_sd;
float	local_variance;

float	variance_explained;
float	variance_explained_array[200];
float	cluster_num_array[200];
float	func_array[200];
float	sd_pointer_array[200];
float	first_diff[200];
float	sec_diff[200];
float	third_diff[200];
float	val_diff;
float	max_third_diff;
float	max_variance_explained;
float	cluster_num_of_max;

int	frames;
int	all_frames;
int	sum_frames = 0;
int	cluster_frames;
int	times;


int matrix_max;
 
FILE *op;
int	cluster_num = 1;

int	cluster_pixels;
float	func;

int	max_cluster_frames = 0;
float	cutoff_number_of_frames = -1;

setlinebuf( stdout );

begin = clock();



/********************************************************************************/
/*										*/
/*Arguments sanity checks, opening file 					*/
/*										*/
/********************************************************************************/



for (i = 0; i < argc; i++)
{
	if ( strncasecmp( argv[i], "-F", 2) == 0 )
	{
		if ( sscanf( argv[i+1], "%d", &FRACTION ) != 1 )
		{
			printf("Error : -fraction expects an integer argument\n");
			printf("Usage : cluster5D [-v] [-fract <integer>] [PCA filename]\n");
			exit(1);
		}
		if ( FRACTION < 0 || FRACTION > 100 )
		{
			printf("Error : argument to -fraction should be an integer between 0 and 100.\n");
			printf("Usage : cluster5D [-v] [-fract <integer>] [PCA filename]\n");
			exit(1);
		}
		AUTO_FRACTION = NO;
		i++;
	}

	if ( strncasecmp (argv[i], "-V", 2) == 0 )
	{
		VERBOSE = YES;
	}
	if ( i == argc - 1)
	{
		if ( i == 0 )
		{
			AUTO_FILENAME = YES;
		}
		else
		{
			fp = fopen(argv[i], "r");
			if (fp == NULL)
 			{
		
  				AUTO_FILENAME = YES;
 			}
		}
	}
}

if ( AUTO_FILENAME == YES )
{
	count = 0;
	count_dpca = 0;
	count_pca = 0;
	if ( (dir = opendir(".")) != NULL )
	{
		while ( (dp = readdir(dir)) != NULL)
		{
			if ( strcmp(dp->d_name,"carma.dPCA.fluctuations.dat") == 0 )
			{
				count_dpca++;
				count++;
			}
			if ( strcmp(dp->d_name,"carma.PCA.fluctuations.dat") == 0 )
			{
				count_pca++;
				count++;
			}
		}
		if ( count == 1 )
		{
			if ( count_dpca == 1)
			{
				fp = fopen("carma.dPCA.fluctuations.dat", "r");
				if (fp == NULL)
				{
					printf ("Error: Cannot open file carma.dPCA.fluctuations.dat\n");
					exit(1);
				}
			}
			if ( count_pca == 1)
			{
				fp = fopen("carma.PCA.fluctuations.dat", "r");
				if (fp == NULL)
				{
					printf ("Error: Cannot open file carma.PCA.fluctuations.dat\n");
					exit(1);
				}
			}	
		}
		else
		if ( count == 2 )
		{
			printf ("Error: There are more than one PCA files in the current folder. Please select one PCA file.\n");
			exit(1); 
		}
		else
		{ 
 			{ 
  				printf ("Error: There is neither an argument for a PCA file nor a PCA file in the current directory.\n");
  				exit(1);
 			}	
		}
	}
}
while ( (fgets(line, sizeof(line), fp)) != NULL )
{
	if ( strlen(line) > 9999)
	{
		printf ("Error: Too big number of columns in the PCA file.\n");
		exit(1);
	}
	if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) != 6 )
	{
			printf ("Error: Invalid file. It must have at least 6 columns.\n");
			exit(1);
	}
}


/************************************************************************************************/
/*												*/
/*First pass to determine limits for each PCA row						*/
/*												*/
/************************************************************************************************/


rewind(fp);
have_limits = 0;
all_frames = 0;
if( VERBOSE == YES)
{
printf("First pass to determine limits ...\n");
printf("Now processing frame ");
}
while( ( fgets (line, sizeof(line), fp) ) != NULL )
{
	if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) == 6 )
	{
		all_frames++;
		if ( VERBOSE == YES)	
		{printf("%8d\b\b\b\b\b\b\b\b",num);}
		if (have_limits == 0)
		{
			pca1_max = pca1;
			pca1_min = pca1;
			pca2_max = pca2;
			pca2_min = pca2;
			pca3_max = pca3;
			pca3_min = pca3;
			pca4_max = pca4;
			pca4_min = pca4;
			pca5_max = pca5;
			pca5_min = pca5;
	
			have_limits = 1;
		}
  
		if ( pca1 > pca1_max )
		pca1_max = pca1;
		if ( pca1 < pca1_min )
		pca1_min = pca1;
		if ( pca2 > pca2_max)
		pca2_max = pca2;
		if ( pca2 < pca2_min)
		pca2_min = pca2;
		if (pca3 > pca3_max)
		pca3_max = pca3;
		if (pca3 < pca3_min)
		pca3_min = pca3;
		if (pca4 > pca4_max)
		pca4_max = pca4;
		if (pca4 <pca4_min)
		pca4_min = pca4;
		if (pca5 > pca5_max)
		pca5_max = pca5;
		if (pca5 < pca5_min)
		pca5_min = pca5;
	}
}
	
if(VERBOSE == YES)
{printf("\n");
printf("%d\tframes will enter the calculation.\n", num);}


/********************************************************************************************************/
/*												        */
/*Second pass to populate the 5-dimensional matrix							*/
/*												        */
/********************************************************************************************************/

DIM = (int) ( pow(all_frames*2, 0.2 ) + 0.5);
rewind (fp);
if(VERBOSE == YES)
{printf("Second pass to populate the 5D matrix ...\n");
printf("Now processing frame ");}
while( ( fgets (line, sizeof(line), fp) ) != NULL )
{
	if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) == 6 )
	{
		if(VERBOSE == YES)
		{printf("%8d\b\b\b\b\b\b\b\b",num);}
		pointer1= (int)( ((pca1 - pca1_min) / (pca1_max - pca1_min)) * DIM + 0.5);
		pointer2= (int)( ((pca2 - pca2_min) / (pca2_max - pca2_min)) * DIM + 0.5);
		pointer3= (int)( ((pca3 - pca3_min) / (pca3_max - pca3_min)) * DIM + 0.5);
		pointer4= (int)( ((pca4 - pca4_min) / (pca4_max - pca4_min)) * DIM + 0.5);
		pointer5= (int)( ((pca5 - pca5_min) / (pca5_max - pca5_min)) * DIM + 0.5);

		original_matrix[pointer1][pointer2][pointer3][pointer4][pointer5]++;
	}
}

if(VERBOSE == YES)
{printf("\n");
fflush(stdout);}

/****************************************/
/*					*/
/*Smoothing					*/
/*					*/
/****************************************/

if(VERBOSE == YES)
{printf("Now smoothing ...");
fflush(stdout);}

for ( i=0 ; i <= DIM ; i++ )
for ( k=0 ; k <= DIM ; k++ )
for ( j=0 ; j <= DIM ; j++ )
for ( l=0 ; l <= DIM ; l++ )
for ( m=0 ; m <= DIM ; m++ )
{
	smoothing(i,k,j,l,m);
}

if(VERBOSE == YES)
{printf("\n");}


/****************************************************************/
/*								*/
/*Calculation of the density threshold				*/
/*								*/
/****************************************************************/







/************************************************************************************************/
/*												*/
/*Calculation of mean, standard deviation and variance						*/
/*												*/
/************************************************************************************************/


N = 0;
sum = 0.0;
mean = 0.0;
for ( i=0 ; i <= DIM ; i++ )
for ( k=0 ; k <= DIM ; k++ )
for ( j=0 ; j <= DIM ; j++ )
for ( l=0 ; l <= DIM ; l++ )
for ( m=0 ; m <= DIM ; m++ )
	{
		value = matrix[i][k][j][l][m];
		if ( value > 0 )
		{
			sum += value;
			N++;
		}
		if ( value > add )
		{
			printf("Error. Increase ... and recompile.\n");
			exit(1);
		}		
	}
mean = sum / N;			


sum_sd = 0.0;
val = 0.0;
for ( i=0 ; i <= DIM ; i++ )
for ( k=0 ; k <= DIM ; k++ )
for ( j=0 ; j <= DIM ; j++ )
for ( l=0 ; l <= DIM ; l++ )
for ( m=0 ; m <= DIM ; m++ )
	{
		value = matrix[i][k][j][l][m];
		if ( value > 0 )
		{
			val = ( value - mean);
			sum_sd += (val * val);
		}
	}
variance = sum_sd / N;
sd = sqrt (variance);

/********************************************************************************************************/
/*													*/
/*Calculation of the density threshold from the given argument						*/
/*													*/
/********************************************************************************************************/

if ( AUTO_FRACTION == NO )
{
	if (VERBOSE == YES)
	{printf ("Testing density threshold: ");
	fflush(stdout);}
	frames =(int) ( (all_frames * FRACTION / 100) + 0.5 );
	frames_count = all_frames;
	times = 0;
	sd_pointer = 0.0;
	while ( frames <= frames_count )
	{
		cutoff = mean + ( sd_pointer * sd );
		frames_count = 0;
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			value1 = original_matrix[i][k][j][l][m];
			value = matrix[i][k][j][l][m];
			if ( (float) value >= cutoff)
			{
				frames_count += value1;
			}
		}
		if (VERBOSE == YES)
		{printf("%10.2f\b\b\b\b\b\b\b\b\b\b", cutoff);
		fflush(stdout);}
		if ( frames >= frames_count )
		{
			break;
		}
		times++;
		sd_pointer = times * 0.1000l;
	}
if (VERBOSE == YES)
{printf("\n");}
if (VERBOSE == YES)
{printf ("Density threshold set to %.2f.\n", cutoff);}
}



/************************************************************************************************************************/
/*															*/
/*Automatically calculation of the density threshold									*/
/*															*/
/************************************************************************************************************************/

if ( AUTO_FRACTION == YES )
{
	for ( i=0 ; i <= DIM ; i++ )
	for ( k=0 ; k <= DIM ; k++ )
	for ( j=0 ; j <= DIM ; j++ )
	for ( l=0 ; l <= DIM ; l++ )
	for ( m=0 ; m <= DIM ; m++ )
	{
		original_matrix[i][k][j][l][m] = matrix[i][k][j][l][m];

	}

	op = fopen ("cluster5D_variance_explained.dat", "w+");
	times = 0;
	sd_pointer = 1.0;
	local_variance = variance ;
	variance_explained = 100 * (local_variance / variance);
	if (VERBOSE == YES)
	{printf ("Testing density threshold: ");
	fflush(stdout);}
	while ( sd_pointer <= 15.0 )
	{
		cutoff = mean + ( sd_pointer * sd );

		if (VERBOSE == YES)
		{printf("%10.2f\b\b\b\b\b\b\b\b\b\b", cutoff);
		fflush(stdout);}
		
		add = MIN_ADD;
		local_sum_sd = 0.0;
		cluster_num = 0;
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			value = matrix[i][k][j][l][m];
			if ( (float) value >= cutoff )
			{
				local_val = (value - mean);
				local_sum_sd += (local_val * local_val);
			}
		}
		
		local_variance = local_sum_sd / N ;		
		
		matrix_max = matrix[0][0][0][0][0];
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			value = matrix[i][k][j][l][m];
			if ( value > matrix_max )
			{
				matrix_max = value;
				pointer1 = i;
				pointer2 = k;
				pointer3 = j;
				pointer4 = l;
				pointer5 = m;
			}
		}
		


		while ( (float) matrix_max > cutoff )
		{
			recursion (pointer1,pointer2,pointer3,pointer4,pointer5);
			cluster_num++;
			matrix_max = matrix[0][0][0][0][0];
			for ( i=0 ; i <= DIM ; i++ )
			for ( k=0 ; k <= DIM ; k++ )
			for ( j=0 ; j <= DIM ; j++ )
			for ( l=0 ; l <= DIM ; l++ )
			for ( m=0 ; m <= DIM ; m++ )
			{
				value = matrix[i][k][j][l][m];
				if ( value > matrix_max && value < MIN_ADD )
				{
					matrix_max = value;
					pointer1 = i;
					pointer2 = k;
					pointer3 = j;
					pointer4 = l;
					pointer5 = m;
				}
			}

			add += 1000000;	
			
		}

		
		variance_explained = 100 * (local_variance / variance);
		func = variance_explained * cluster_num; 
		func_array[times] = func;
		sd_pointer_array[times] = sd_pointer;
		variance_explained_array[times] = variance_explained;
		cluster_num_array[times] = cluster_num;
		fprintf(op,"%4d\t%7.4f\t%7.2f\n",cluster_num,variance_explained,cutoff);
		times++;
		sd_pointer += 0.1000l;
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			matrix[i][k][j][l][m] = original_matrix[i][k][j][l][m];

		}
	}


	for ( i = 0 ; i <= times - 1 ; i++)
	{
		val_diff = ( func_array[i+1] - func_array[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
		first_diff[i] = val_diff;
	}

	for ( i = 0 ; i <= times - 2 ; i++)
	{
		val_diff = ( first_diff[i+1] - first_diff[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
		sec_diff[i] = val_diff;
	}

	for ( i = 0 ; i <= times - 3 ; i++)
	{
		val_diff = ( sec_diff[i+1] - sec_diff[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
		third_diff[i] = val_diff;
	}

	max_third_diff = third_diff[0];
	sd_pointer = sd_pointer_array[0];
	max_variance_explained = 0.0;
	cluster_num_of_max = 0;	
	
	for ( i = 0 ; i <= times - 3 ; i++)
	{	
		val_diff = third_diff[i];
		if ( val_diff > max_third_diff )
		{
			max_third_diff = val_diff;
			cluster_num_of_max = cluster_num_array[i];
		}
	}
			
	for ( j = 0 ; j <= times - 3 ; j++)
	{
		if ( ( cluster_num_array[j] == cluster_num_of_max ) && ( variance_explained_array[j] > max_variance_explained ) )
				{
					max_variance_explained = variance_explained_array[j];
					
					sd_pointer = sd_pointer_array[j];
				}
	}
	
	cutoff =  mean + ( sd_pointer * sd );
	if (VERBOSE == YES)
	{printf("\n");}
	if (VERBOSE == YES)
	{printf ("Density threshold set to %.2f.\n", cutoff);}

	fclose(op);

	for ( i=0 ; i <= DIM ; i++ )
	for ( k=0 ; k <= DIM ; k++ )
	for ( j=0 ; j <= DIM ; j++ )
	for ( l=0 ; l <= DIM ; l++ )
	for ( m=0 ; m <= DIM ; m++ )
	{
		matrix[i][k][j][l][m] = original_matrix[i][k][j][l][m];

	}

}


cluster_num = 1;

/************************/
/*			*/			
/*Clustering		*/
/*			*/
/************************/



/****************************************************************/
/*								*/
/*Finding the pixel with the maximum value			*/
/*								*/
/****************************************************************/


matrix_max = matrix[0][0][0][0][0];
	for ( i=0 ; i <= DIM ; i++ )
	for ( k=0 ; k <= DIM ; k++ )
	for ( j=0 ; j <= DIM ; j++ )
	for ( l=0 ; l <= DIM ; l++ )
	for ( m=0 ; m <= DIM ; m++ )
	{
		value = matrix[i][k][j][l][m];
		if ( value > matrix_max )
		{
			matrix_max = value;
			pointer1 = i;
			pointer2 = k;
			pointer3 = j;
			pointer4 = l;
			pointer5 = m;
		}
	}



/************************************************************************/
/*									*/
/*Creating the output files and clustering using the recursive function	*/
/*									*/
/************************************************************************/

op = fopen ("carma.5-D.clusters.dat", "w+");
if (VERBOSE == YES)
{printf("Clustering now ...\n");}

while ( (float) matrix_max > cutoff )
{
	
	cluster_pixels = 0;
	cluster_frames = 0;
	recursion (pointer1,pointer2,pointer3,pointer4,pointer5);


/****************************************************************************************/
/*											*/
/* PCA file pass to write the frames that belong to this cluster			*/
/*and to calculate the percentage of pixels and frames					*/
/*											*/
/****************************************************************************************/



	rewind(fp);
	while( ( fgets (line, sizeof(line), fp) ) != NULL )
	{
		if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) == 6 )
		{
			pointer1= (int)( ((pca1 - pca1_min) / (pca1_max - pca1_min)) * DIM + 0.5);
			pointer2= (int)( ((pca2 - pca2_min) / (pca2_max - pca2_min)) * DIM + 0.5);
			pointer3= (int)( ((pca3 - pca3_min) / (pca3_max - pca3_min)) * DIM + 0.5);
			pointer4= (int)( ((pca4 - pca4_min) / (pca4_max - pca4_min)) * DIM + 0.5);
			pointer5= (int)( ((pca5 - pca5_min) / (pca5_max - pca5_min)) * DIM + 0.5);
	
			if ( -matrix[pointer1][pointer2][pointer3][pointer4][pointer5] >= add + cutoff )
			{
				fprintf (op,"%10d  %10d  %10f   %10f   %10f   %10f   %10f\n", num,cluster_num,pca1,pca2,pca3,pca4,pca5);
				cluster_frames++;
			}
	 	}
	}
	for ( i=0 ; i <= DIM ; i++ )
	for ( k=0 ; k <= DIM ; k++ )
	for ( j=0 ; j <= DIM ; j++ )
	for ( l=0 ; l <= DIM ; l++ )
	for ( m=0 ; m <= DIM ; m++ )
	{
		value = -matrix[i][k][j][l][m];
		if ( value >= add + cutoff )
		{
			cluster_pixels++;
		}
	}
	if ( cluster_frames != 0)
	{
		if (VERBOSE == YES)
		{printf ("Cluster %5d located, contains %10d frames.\n", cluster_num, cluster_frames);}
		sum_frames += cluster_frames;
		if ( cluster_frames > max_cluster_frames )
		{
			max_cluster_frames = cluster_frames;
		}
		if ( cluster_num == 20 )
		{
			cutoff_number_of_frames = max_cluster_frames / 10000;
		}
		if ( cutoff_number_of_frames > 0 && cluster_frames < cutoff_number_of_frames )
		{
			if (VERBOSE == YES)
			{printf ("With several small (<%d frames) clusters following ...\n", (int) cutoff_number_of_frames);}
			end = clock();
			time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
			if (VERBOSE == YES)
			{printf("All done in %.1f minutes.\n", time_spent / 60);}
			exit(1);
		}		

		cluster_num++;
	}	

	
	


/****************************************************************************************/
/*											*/				
/*Finding the pixel with the next maximum value	and repeating				*/
/*the above steps until the value of pixel reaches the threshold			*/
/*											*/
/****************************************************************************************/


	matrix_max = matrix[0][0][0][0][0];
	for ( i=0 ; i <= DIM ; i++ )
	for ( k=0 ; k <= DIM ; k++ )
	for ( j=0 ; j <= DIM ; j++ )
	for ( l=0 ; l <= DIM ; l++ )
	for ( m=0 ; m <= DIM ; m++ )
	{
		value = matrix[i][k][j][l][m];
		if ( value > matrix_max && value < MIN_ADD )
			{
				matrix_max = value;
				pointer1 = i;
				pointer2 = k;
				pointer3 = j;
				pointer4 = l;
				pointer5 = m;
			}
	}
	add += 1000000;
}
	fclose(op);
	if ( AUTO_FRACTION == YES )
	{
		if (VERBOSE == YES)
		{printf("%.2f%% of frames have been clustered.\n", 100.0 * sum_frames / all_frames);}
	}	










	/************************************************************************/
	/*									*/
	/*If the percentage of frames assigned to clusters is less than 10%,	*/
	/*repeat the whole procedure without smoothing				*/ 
	/*									*/
	/************************************************************************/
	
	if ( (100.0 * sum_frames / all_frames) < 10 )
	{
		if (VERBOSE == YES) 
		{printf("Too few frames assigned to clusters.\nRepeating the procedure without smoothing:\n");}
		
		sum_frames = 0;		
	
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			original_matrix[i][k][j][l][m] = 0;
		}

		/********************************************************************************************************/
		/*												        */
		/*Second pass to populate the 5-dimensional matrix							*/
		/*												        */
		/********************************************************************************************************/

		DIM = (int) ( pow(all_frames*2, 0.2 ) + 0.5);
		rewind (fp);
		if(VERBOSE == YES)
		{printf("Second pass to populate the 5D matrix ...\n");
		printf("Now processing frame ");}
		while( ( fgets (line, sizeof(line), fp) ) != NULL )
		{
			if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) == 6 )
			{
				if(VERBOSE == YES)
				{printf("%8d\b\b\b\b\b\b\b\b",num);}
				pointer1= (int)( ((pca1 - pca1_min) / (pca1_max - pca1_min)) * DIM + 0.5);
				pointer2= (int)( ((pca2 - pca2_min) / (pca2_max - pca2_min)) * DIM + 0.5);
				pointer3= (int)( ((pca3 - pca3_min) / (pca3_max - pca3_min)) * DIM + 0.5);
				pointer4= (int)( ((pca4 - pca4_min) / (pca4_max - pca4_min)) * DIM + 0.5);
				pointer5= (int)( ((pca5 - pca5_min) / (pca5_max - pca5_min)) * DIM + 0.5);

				original_matrix[pointer1][pointer2][pointer3][pointer4][pointer5]++;
			}
		}

		if(VERBOSE == YES)
		{printf("\n");
		fflush(stdout);}

		
		
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			matrix[i][k][j][l][m] = original_matrix[i][k][j][l][m];
		}


		/****************************************************************/
		/*								*/
		/*Calculation of the density threshold				*/
		/*								*/
		/****************************************************************/

	



		/************************************************************************************************/
		/*												*/
		/*Calculation of mean, standard deviation and variance						*/
		/*												*/
		/************************************************************************************************/




		N = 0;
		sum = 0.0;
		mean = 0.0;
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			value = matrix[i][k][j][l][m];
			if ( value > 0 )
			{
				sum += value;
				N++;
			}
			if ( value > add )
			{
				printf("Error. Increase ... and recompile.\n");
				exit(1);
			}		
		}
		mean = sum / N;			


		sum_sd = 0.0;
		val = 0.0;
		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			value = matrix[i][k][j][l][m];
			if ( value > 0 )
			{
				val = ( value - mean);
				sum_sd += (val * val);
			}
		}
		variance = sum_sd / N;
		sd = sqrt (variance);

		/********************************************************************************************************/
		/*													*/
		/*Calculation of the density threshold from the given argument						*/
		/*													*/
		/********************************************************************************************************/

		if ( AUTO_FRACTION == NO )
		{
			if (VERBOSE == YES)
			{printf ("Testing density threshold: ");
			fflush(stdout);}
			frames =(int) ( (all_frames * FRACTION / 100) + 0.5 );
			frames_count = all_frames;
			times = 0;
			sd_pointer = 0.0;
			while ( frames <= frames_count )
			{
				cutoff = mean + ( sd_pointer * sd );
				frames_count = 0;
				for ( i=0 ; i <= DIM ; i++ )
				for ( k=0 ; k <= DIM ; k++ )
				for ( j=0 ; j <= DIM ; j++ )
				for ( l=0 ; l <= DIM ; l++ )
				for ( m=0 ; m <= DIM ; m++ )
				{
					value1 = original_matrix[i][k][j][l][m];
					value = matrix[i][k][j][l][m];
					if ( (float) value >= cutoff)
					{
						frames_count += value1;
					}
				}
				if (VERBOSE == YES)
				{printf("%10.2f\b\b\b\b\b\b\b\b\b\b", cutoff);
				fflush(stdout);}
				if ( frames >= frames_count )
				{
					break;
				}
				times++;
				sd_pointer = times * 0.1000l;
			}
		if (VERBOSE == YES)
		{printf("\n");}
		if (VERBOSE == YES)
		{printf ("Density threshold set to %.2f.\n", cutoff);}
		}



		/****************************************************************************************************************/
		/*														*/
		/*Automatically calculation of the density threshold								*/
		/*														*/
		/****************************************************************************************************************/

		if ( AUTO_FRACTION == YES )
		{
			op = fopen ("cluster5D_variance_explained.dat", "w+");
			times = 0;
			sd_pointer = 0.0;
			local_variance = variance ;
			variance_explained = 100 * (local_variance / variance);
			if (VERBOSE == YES)
			{printf ("Testing density threshold: ");
			fflush(stdout);}
			while ( sd_pointer <= 15.0 )
			{
				cutoff = mean + ( sd_pointer * sd );

				if (VERBOSE == YES)
				{printf("%10.2f\b\b\b\b\b\b\b\b\b\b", cutoff);
				fflush(stdout);}
		
				add = MIN_ADD;
				local_sum_sd = 0.0;
				cluster_num = 0;
				for ( i=0 ; i <= DIM ; i++ )
				for ( k=0 ; k <= DIM ; k++ )
				for ( j=0 ; j <= DIM ; j++ )
				for ( l=0 ; l <= DIM ; l++ )
				for ( m=0 ; m <= DIM ; m++ )
				{
					value = matrix[i][k][j][l][m];
					if ( (float) value >= cutoff )
					{
						local_val = (value - mean);
						local_sum_sd += (local_val * local_val);
					}
				}
		
				local_variance = local_sum_sd / N ;		
		
				matrix_max = matrix[0][0][0][0][0];
				for ( i=0 ; i <= DIM ; i++ )
				for ( k=0 ; k <= DIM ; k++ )
				for ( j=0 ; j <= DIM ; j++ )
				for ( l=0 ; l <= DIM ; l++ )
				for ( m=0 ; m <= DIM ; m++ )
				{
					value = matrix[i][k][j][l][m];
					if ( value > matrix_max )
					{
						matrix_max = value;
						pointer1 = i;
						pointer2 = k;
						pointer3 = j;
						pointer4 = l;
						pointer5 = m;
					}
				}
		


				while ( (float) matrix_max > cutoff )
				{
					recursion (pointer1,pointer2,pointer3,pointer4,pointer5);
					cluster_num++;
					matrix_max = matrix[0][0][0][0][0];
					for ( i=0 ; i <= DIM ; i++ )
					for ( k=0 ; k <= DIM ; k++ )
					for ( j=0 ; j <= DIM ; j++ )
					for ( l=0 ; l <= DIM ; l++ )
					for ( m=0 ; m <= DIM ; m++ )
					{
						value = matrix[i][k][j][l][m];
						if ( value > matrix_max && value < MIN_ADD )
						{
							matrix_max = value;
							pointer1 = i;
							pointer2 = k;
							pointer3 = j;
							pointer4 = l;
							pointer5 = m;
						}
					}
		
					add += 1000000;	
					
				}

		
				variance_explained = 100 * (local_variance / variance);
				func = variance_explained * cluster_num; 
				func_array[times] = func;
				sd_pointer_array[times] = sd_pointer;
				variance_explained_array[times] = variance_explained;
				cluster_num_array[times] = cluster_num;
				fprintf(op,"%4d\t%7.4f\t%7.2f\n",cluster_num,variance_explained,cutoff);
				times++;
				sd_pointer += 0.1000l;
				for ( i=0 ; i <= DIM ; i++ )
				for ( k=0 ; k <= DIM ; k++ )
				for ( j=0 ; j <= DIM ; j++ )
				for ( l=0 ; l <= DIM ; l++ )
				for ( m=0 ; m <= DIM ; m++ )
				{
					matrix[i][k][j][l][m] = original_matrix[i][k][j][l][m];
		
				}
			}		


			for ( i = 0 ; i <= times - 1 ; i++)
			{
				val_diff = ( func_array[i+1] - func_array[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
				first_diff[i] = val_diff;
			}

			for ( i = 0 ; i <= times - 2 ; i++)
			{
				val_diff = ( first_diff[i+1] - first_diff[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
				sec_diff[i] = val_diff;
			}

			for ( i = 0 ; i <= times - 3 ; i++)
			{
				val_diff = ( sec_diff[i+1] - sec_diff[i] ) / ( sd_pointer_array[i+1] - sd_pointer_array[i] );
				third_diff[i] = val_diff;
			}

			max_third_diff = third_diff[0];
			sd_pointer = sd_pointer_array[0];
			max_variance_explained = 0.0;
			cluster_num_of_max = 0;	
	
			for ( i = 0 ; i <= times - 3 ; i++)
			{	
				val_diff = third_diff[i];
				if ( val_diff > max_third_diff )
				{
					max_third_diff = val_diff;
					cluster_num_of_max = cluster_num_array[i];
				}
			}
			
			for ( j = 0 ; j <= times - 3 ; j++)
			{
				if ( ( cluster_num_array[j] == cluster_num_of_max ) && ( variance_explained_array[j] > max_variance_explained ) )
				{
					max_variance_explained = variance_explained_array[j];
					
					sd_pointer = sd_pointer_array[j];
				}
			}
	
			cutoff =  mean + ( sd_pointer * sd );
			if (VERBOSE == YES)
			{printf("\n");}
			if (VERBOSE == YES)
			{printf ("Density threshold set to %.2f.\n", cutoff);}

			fclose(op);
		}

		for ( i=0 ; i <= DIM ; i++ )
		for ( k=0 ; k <= DIM ; k++ )
		for ( j=0 ; j <= DIM ; j++ )
		for ( l=0 ; l <= DIM ; l++ )
		for ( m=0 ; m <= DIM ; m++ )
		{
			matrix[i][k][j][l][m] = original_matrix[i][k][j][l][m];
		
		}

		cluster_num = 1;


		/************************/
		/*			*/			
		/*Clustering		*/
		/*			*/
		/************************/



		/****************************************************************/
		/*								*/
		/*Finding the pixel with the maximum value			*/
		/*								*/
		/****************************************************************/


		matrix_max = matrix[0][0][0][0][0];
			for ( i=0 ; i <= DIM ; i++ )
			for ( k=0 ; k <= DIM ; k++ )
			for ( j=0 ; j <= DIM ; j++ )
			for ( l=0 ; l <= DIM ; l++ )
			for ( m=0 ; m <= DIM ; m++ )
			{
				value = matrix[i][k][j][l][m];
				if ( value > matrix_max )
				{
					matrix_max = value;
					pointer1 = i;
					pointer2 = k;
					pointer3 = j;
					pointer4 = l;
					pointer5 = m;
				}
			}
	


		/************************************************************************/
		/*									*/
		/*Creating the output files and clustering using the recursive function	*/
		/*									*/
		/************************************************************************/

		op = fopen ("carma.5-D.clusters.dat", "w+");
		if (VERBOSE == YES)
		{printf("Clustering now ...\n");}

		while ( (float) matrix_max > cutoff )
		{
	
			cluster_pixels = 0;
			cluster_frames = 0;
			recursion (pointer1,pointer2,pointer3,pointer4,pointer5);
		

		/****************************************************************************************/
		/*											*/
		/* PCA file pass to write the frames that belong to this cluster			*/
		/*and to calculate the percentage of pixels and frames					*/
		/*											*/
		/****************************************************************************************/



			rewind(fp);
			while( ( fgets (line, sizeof(line), fp) ) != NULL )
			{
				if ( (sscanf(line,"%d   %f   %f   %f   %f   %f",&num,&pca1,&pca2,&pca3,&pca4,&pca5)) == 6 )
				{
					pointer1= (int)( ((pca1 - pca1_min) / (pca1_max - pca1_min)) * DIM + 0.5);
					pointer2= (int)( ((pca2 - pca2_min) / (pca2_max - pca2_min)) * DIM + 0.5);
					pointer3= (int)( ((pca3 - pca3_min) / (pca3_max - pca3_min)) * DIM + 0.5);
					pointer4= (int)( ((pca4 - pca4_min) / (pca4_max - pca4_min)) * DIM + 0.5);
					pointer5= (int)( ((pca5 - pca5_min) / (pca5_max - pca5_min)) * DIM + 0.5);
		
					if ( -matrix[pointer1][pointer2][pointer3][pointer4][pointer5] >= add + cutoff )
					{
						fprintf (op,"%10d  %10d  %10f   %10f   %10f   %10f   %10f\n", num,cluster_num,pca1,pca2,pca3,pca4,pca5);
						cluster_frames++;
					}
			 	}
			}
			for ( i=0 ; i <= DIM ; i++ )
			for ( k=0 ; k <= DIM ; k++ )
			for ( j=0 ; j <= DIM ; j++ )
			for ( l=0 ; l <= DIM ; l++ )
			for ( m=0 ; m <= DIM ; m++ )
			{
				value = -matrix[i][k][j][l][m];
				if ( value >= add + cutoff )
				{
					cluster_pixels++;
				}
			}
	
			if (VERBOSE == YES)
			{printf ("Cluster %5d located, contains %10d frames.\n", cluster_num, cluster_frames);}
			sum_frames += cluster_frames;
			if ( cluster_frames > max_cluster_frames )
			{
				max_cluster_frames = cluster_frames;
			}
			if ( cluster_num == 20 )
			{
				cutoff_number_of_frames = max_cluster_frames / 10000;
			}
			if ( cutoff_number_of_frames > 0 && cluster_frames < cutoff_number_of_frames )
			{
				if (VERBOSE == YES)
				{printf ("With several small (<%d frames) clusters following ...\n", (int) cutoff_number_of_frames);}
				end = clock();
				time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
				if (VERBOSE == YES)
				{printf("All done in %.1f minutes.\n", time_spent / 60);}
				exit(1);
			}
			cluster_num++;


		/****************************************************************************************/
		/*										     */				
		/*Finding the pixel with the next maximum value	and repeating				*/
		/*the above steps until the value of pixel reaches the threshold			*/
		/*											*/
		/****************************************************************************************/


			matrix_max = matrix[0][0][0][0][0];
			for ( i=0 ; i <= DIM ; i++ )
			for ( k=0 ; k <= DIM ; k++ )
			for ( j=0 ; j <= DIM ; j++ )
			for ( l=0 ; l <= DIM ; l++ )
			for ( m=0 ; m <= DIM ; m++ )
			{
				value = matrix[i][k][j][l][m];
				if ( value > matrix_max && value < MIN_ADD )
				{
					matrix_max = value;
					pointer1 = i;
					pointer2 = k;
					pointer3 = j;
					pointer4 = l;
					pointer5 = m;
				}
			}
			add += 1000000;
		}
		fclose(op);
		if ( AUTO_FRACTION == YES )
		{
			if (VERBOSE == YES)
			{printf("%.2f%% of frames have been clustered.\n", 100.0 * sum_frames / all_frames);}
		}	






	}



	/********************************/
	/*				*/
	/*End of repeating procedure	*/
	/*				*/
	/********************************/




	




	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
		
	if (VERBOSE == YES)
	{printf("All done in %.1f minutes.\n", time_spent / 60);}
	return 0;
}


/****************************************************************************************/
/*											*/
/*The recursive function "recursion" checks the value of pixels that			*/
/*are around the pixel with the maximum value. For these pixels that have		*/
/*value >= threshold, it checks recursively the value of their around pixels		*/
/*											*/
/****************************************************************************************/




void recursion(int p1, int p2, int p3, int p4, int p5)
	{
		int i,k,j,l,m;
		int	val;
		for ( i = p1 - 1; i <= p1 + 1; i++ )
		for ( k = p2 - 1; k <= p2 + 1; k++ )
		for ( j = p3 - 1; j <= p3 + 1; j++ )
		for ( l = p4 - 1; l <= p4 + 1; l++ )
		for ( m = p5 - 1; m <= p5 + 1; m++ )
		{
		if( i >= 0 && i <= DIM && k >= 0 && k <= DIM && j >= 0 && j <= DIM && l >= 0 && l <= DIM && m >= 0 && m <= DIM )
		{
			val = matrix[i][k][j][l][m];
			if ( (float) val >= cutoff && val < add )
			{
				val += add;
				matrix[i][k][j][l][m] = -val;
				recursion(i,k,j,l,m);
			}
		}
		}
	}

/****************************************************************/
/*								*/
/*Smoothing function						*/
/*								*/
/****************************************************************/

void smoothing(int p1, int p2, int p3, int p4, int p5)
{
	int i,k,j,l,m;
	int val;
	int sum, count, average;
	sum = 0;
	count = 0;
	for ( i = p1 - 1; i <= p1 + 1; i++ )
	for ( k = p2 - 1; k <= p2 + 1; k++ )
	for ( j = p3 - 1; j <= p3 + 1; j++ )
	for ( l = p4 - 1; l <= p4 + 1; l++ )
	for ( m = p5 - 1; m <= p5 + 1; m++ )
	{
		if( i >= 0 && i <= DIM && k >= 0 && k <= DIM && j >= 0 && j <= DIM && l >= 0 && l <= DIM && m >= 0 && m <= DIM )
		{
			val = original_matrix[i][k][j][l][m];
			sum += val;
			count++;
		}
	}
	average = (int) ( sum / count + 0.5);
	matrix[p1][p2][p3][p4][p5] = average;
}
