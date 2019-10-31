#include "predict_header.h"
using namespace Rcpp;

// [[Rcpp::export]]
void ReleaseExondata(Exondata target)
{
	std::vector<std::string>().swap(target.sample_name);
	std::vector<std::string>().swap(target.TC_name);
	std::vector<std::vector<double> >().swap(target.data);
}

// [[Rcpp::export]]
//check the input data format
int CheckTCid(char **lib_TC, std::vector<std::string> in_TC, int Length)  
{
	for (int i = 0; i < Length; i++)
	{
		if (strcmp(lib_TC[i], in_TC[i].c_str()) != 0)
		{
			return 1;
		}
	}

	return 0;
}

// [[Rcpp::export]]
//match input data with separate character
int MatchExon_sep(char **lib_TC, Exondata *indata, int Length, int *match_idx, std::string sep_str)
{
    int error_flag = 1;
    std::string id_str;
    std::size_t sep_pos;
    std::string temp_id;
    std::unordered_map<std::string, int> id_map;
    
    for (int i = 0; i < indata->TC_name.size(); i++)
    {
        id_str = indata->TC_name[i].c_str();
        sep_pos = id_str.find(sep_str);
        temp_id = id_str.substr(0, sep_pos);
        id_map[temp_id] = i+1;
    }
    
    for (int i = 0; i < Length; i++)
    {
        id_str = lib_TC[i];
        sep_pos = id_str.find(sep_str);
        temp_id = id_str.substr(0, sep_pos);
        
        match_idx[i] = id_map[temp_id] - 1;
    }
    
    for (int i = 0; i < Length; i++)
    {
        if(match_idx[i] != -1)
        {
            error_flag = 0;
            break;
        }
    }
    
    if(error_flag == 1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// [[Rcpp::export]]
//match data with extract match
int MatchExon(char **lib_TC, Exondata *indata, int Length, int *match_idx)
{
    int error_flag = 1;
    std::unordered_map<std::string, int> id_map;
    
    for (int i = 0; i < indata->TC_name.size(); i++)
    {
        id_map[indata->TC_name[i].c_str()] = i+1;
    }
    
    for (int i = 0; i < Length; i++)
    {
        match_idx[i] = id_map[lib_TC[i]] - 1;
    }
    
    for (int i = 0; i < Length; i++)
    {
        if(match_idx[i] != -1)
        {
            error_flag = 0;
            break;
        }
    }
    
    if(error_flag == 1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// [[Rcpp::export]]
//write matched gene expression data matrix
int WriteExpr(double **data_out, std::vector<std::string> outname, char **lib_TC, char *outfile, int Length, int sample_size)
{
    FILE * pFile;
    
    pFile = fopen(outfile, "w");
    if (pFile != NULL)
    {
        fprintf(pFile, "%s\t", "gene_id");
        
        for (int i = 0; i < sample_size-1; i++)
        {
            fprintf(pFile, "%s\t", outname[i].c_str());
        }
        
        fprintf(pFile, "%s\n", outname[sample_size-1].c_str());
        
        for (int i = 0; i < Length; i++)
        {
            fprintf(pFile, "%s\t", lib_TC[i]);
            
            for (int j = 0; j < sample_size-1; j++)
            {
                fprintf(pFile, "%lf\t", pow(2, data_out[j][i]) - 1);
            }
            
            fprintf(pFile, "%lf\n", pow(2, data_out[sample_size-1][i]) - 1);
            
        }
        fclose(pFile);
    }
    else
    {
        std::cout << "Error! Cannot write file " << outfile << std::endl;
        return 1;
    }

    return 0;
}

// [[Rcpp::export]]
//read in gene expression data
int ReadinExon(char filename[255], Exondata *indata)  
{
	std::ifstream infile(filename);
	std::string line_temp;
	std::string line_token;

	if (infile.is_open())
	{
		std::cout << "Reading input file " << filename << std::endl;
		std::getline(infile, line_temp);
                if(line_temp[line_temp.length()-1]=='\r') {line_temp[line_temp.length()-1] = '\0';}     //handle windows file
		std::stringstream each_line(line_temp);
		std::getline(each_line, line_token, '\t');          //ignore the first string 'TranscriptCluster'
		while (std::getline(each_line, line_token, '\t'))
		{
			if (line_token.empty()){ break;}
			indata->sample_name.push_back(line_token);
		}

		while (!infile.eof())
		{
			std::getline(infile, line_temp);
                        if(line_temp[line_temp.length()-1]=='\r') {line_temp[line_temp.length()-1] = '\0';}
			if (line_temp.empty()){ break; }
			std::stringstream each_line(line_temp);
			std::getline(each_line, line_token, '\t');
			indata->TC_name.push_back(line_token);
			std::vector<double> data_temp;
			while (std::getline(each_line, line_token, '\t'))
			{
				if (line_token.empty()){ break;}
				data_temp.push_back(atof(line_token.c_str()));
			}
			indata->data.push_back(data_temp);
		}

		infile.close();
	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}

	return 0;

}

// [[Rcpp::export]]
//shell sorting algorithm
void ShellSort(double *num, int *index, int numLength)    
{
	int i, j, increment, temp_idx;
	double temp;

	for (increment = numLength / 2; increment > 0; increment /= 2)
	{
		for (i = increment; i<numLength; i++)
		{
			temp = num[i];
			temp_idx = index[i];
			for (j = i; j >= increment; j -= increment)
			{

				if (temp > num[j - increment])
				{
					num[j] = num[j - increment];
					index[j] = index[j - increment];
				}
				else
				{
					break;
				}
			}
			num[j] = temp;
			index[j] = temp_idx;
		}
	}
	return;
}
// [[Rcpp::export]]
//get ranking of the sorted data
void GetRank(double *rank, double *num,int numLength)
{
    int i,j,k;
    i = 0;
    while (i < numLength)
    {
        j = i;
        while ((j < numLength - 1) && (num[j]  == num[j + 1]))
        {
            j++;
        }
        if (i != j)
        {
            for (k = i; k <= j; k++)
            {
                rank[k] = (i + j + 2) / 2.0;
            }
        }
        else
        {
            rank[i] = i + 1;
        }
        i = j + 1;
    }
    return;
}

// [[Rcpp::export]]
//quantile normalizaiton with a vector of known quantile
void QuantileNorm(double *indata, double *quantile, int dataLength) 
{
	int *index = new int[dataLength];
	double *indata_copy = new double[dataLength];
    double *ranks = new double[dataLength];

	for (int i = 0; i < dataLength; i++)
	{
		index[i] = i;
		indata_copy[i] = indata[i];
	}

	ShellSort(indata_copy, index, dataLength);
    GetRank(ranks, indata_copy, dataLength);
    
	for (int i = 0; i < dataLength; i++)
	{
        if (ranks[i] - floor(ranks[i]) > 0.4)
        {
            indata[index[i]] = 0.5 * (quantile[(size_t)floor(ranks[i])-1] + quantile[(size_t)floor(ranks[i])]);
        }
        else
        {
            indata[index[i]] = quantile[(size_t)floor(ranks[i])-1];
        }
	}

	delete[] index;
    delete[] ranks;
	delete[] indata_copy;
	return;
}

// [[Rcpp::export]]
//stadardization
void StandardizeRow(double **data_in, double *Mean, double *SD, int Length, int sample_size)  
{
	for (int i = 0; i < Length; i++)
	{
		if (SD[i] != 0)
		{
			for (int j = 0; j < sample_size; j++)
			{
				data_in[j][i] = (data_in[j][i] - Mean[i]) / SD[i];
			}
		}
		else
		{
			for (int j = 0; j < sample_size; j++)
			{
				data_in[j][i] = 0;
			}
		}

	}
	return;
}

// [[Rcpp::export]]
//Reverse standardization
void StandardizeRow_r(double **data_in, double *Mean, double *SD, int Length, int sample_size)  
{
	for (int i = 0; i < Length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			data_in[j][i] = data_in[j][i] * SD[i] + Mean[i];
		}

	}
	return;
}

// [[Rcpp::export]]
//calculate average value within each gene cluster
void ClusterMean(double **data_matrix, double **data_mean, int *cluster_idx, int p_length, int c_length, int sample_size)
{
	double *data_sum = new double[sample_size];
	int data_count;

	for (int i = 0; i < c_length; i++)
	{
		for (int m = 0; m < sample_size; m++)
		{
			data_sum[m] = 0;
		}
		data_count = 0;
		for (int j = 0; j < p_length; j++)
		{
			if (cluster_idx[j] == i + 1)                      
			{
				for (int k = 0; k < sample_size; k++)
				{
					data_sum[k] += data_matrix[k][j];
				}
				data_count++;
			}
		}
		for (int m = 0; m < sample_size; m++)
		{
			data_mean[m][i] = data_sum[m] / (double)data_count;
		}

	}

	delete[] data_sum;
	return;
}

// [[Rcpp::export]]
//read in the pre-built prediction model file
int ReadinModel(char filename[255], double *quantile_in, double *exon_mean, double *exon_sd, double **coef, double *DNase_mean, double *DNase_sd, int **pre_idx, char **TC_id, int *cluster_idx, char **select_loci, int p_length, int var_length, int loci_length, double **dis_matrix, int **DH_cluster, double **DH_coef1, double **DH_coef2, double **DH_coef3, int **DH_pre_idx1, int **DH_pre_idx2, int **DH_pre_idx3, int DH_num1, int DH_num2, int DH_num3)
{
	int dump[8];
	FILE *pFile;
	pFile = fopen(filename, "rb");
	if (pFile != NULL)
	{
		std::cout << "Reading library file " << filename << std::endl;
		fread(dump, sizeof(int), 8, pFile);
		fread(quantile_in, sizeof(double), p_length, pFile);
		fread(exon_mean, sizeof(double), p_length, pFile);
		fread(exon_sd, sizeof(double), p_length, pFile);
		for (int k = 0; k < var_length; k++)
		{
			fread(coef[k], sizeof(double), loci_length, pFile);
		}
		fread(DNase_mean, sizeof(double), loci_length, pFile);
		fread(DNase_sd, sizeof(double), loci_length, pFile);
		for (int k = 0; k < var_length; k++)
		{
			fread(pre_idx[k], sizeof(int), loci_length, pFile);
		}
		fread(cluster_idx, sizeof(int), p_length, pFile);
		for (int k = 0; k < loci_length; k++)
		{
			fread(select_loci[k], sizeof(char), 30, pFile);
		}
		for (int k = 0; k < p_length; k++)
		{
			fread(TC_id[k], sizeof(char), 30, pFile);
		}
		for (int k = 0; k < 3; k++)
		{
			fread(dis_matrix[k], sizeof(double), loci_length, pFile);
		}
		for (int k = 0; k < 3; k++)
		{
			fread(DH_cluster[k], sizeof(int), loci_length, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef1[k], sizeof(double), DH_num1, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef2[k], sizeof(double), DH_num2, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef3[k], sizeof(double), DH_num3, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx1[k], sizeof(int), DH_num1, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx2[k], sizeof(int), DH_num2, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx3[k], sizeof(int), DH_num3, pFile);
		}
	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}

	fclose(pFile);
	return 0;
}

// [[Rcpp::export]]
//read the parameters used in the prediction model
int ReadPar(char filename[255], int &loci_size, int &predictor_size, int &cluster_size, int &bin_size, int &var_size, int &DH_num1, int &DH_num2, int &DH_num3)
{
	FILE *pFile;
	int param[8];
	pFile = fopen(filename, "rb");
	if (pFile != NULL)
	{
		std::cout << "Reading model parameter..." << std::endl;
		fread(param, sizeof(int), 8, pFile);
		loci_size = param[0];
		predictor_size = param[1];
		cluster_size = param[2];
		bin_size = param[3];
		var_size = param[4];
		DH_num1 = param[5];
		DH_num2 = param[6];
		DH_num3 = param[7];
		fclose(pFile);

	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}
	return 0;
}

// [[Rcpp::export]]
//regression according to known coefficients and predictor indexes
void Regression(double **predictor, double **output, double **coef, int **predictor_idx, int var_length, int loci_length, int sample_size)
{
	for (int i = 0; i < loci_length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			output[j][i] = 0;
			for (int k = 0; k < var_length; k++)
			{
				output[j][i] += coef[k][i] * predictor[j][predictor_idx[k][i] - 1];               //index from 1 in library file.
			}

		}
	}
	return;
}

// [[Rcpp::export]]
//model average from different level of prediction
void ModelAverage(double **output, double **DH_pre1, double **DH_pre2, double **DH_pre3, double **dis_matrix, int **DH_cluster, int loci_length, int sample_size)
{
	double weight;
	for (int i = 0; i < loci_length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			weight = dis_matrix[0][i] + dis_matrix[1][i] + dis_matrix[2][i] + 1;
			output[j][i] = (output[j][i] + DH_pre1[j][DH_cluster[0][i] - 1] * dis_matrix[0][i] + DH_pre2[j][DH_cluster[1][i] - 1] * dis_matrix[1][i] + DH_pre3[j][DH_cluster[2][i] - 1] * dis_matrix[2][i]) / weight ;               //index from 1 in library file.
		}
	}
	return;
}

// [[Rcpp::export]]
//write ouput file in the format of data matrix or wig file
int WriteWIG(double **data_out, char **select_idx, std::vector<std::string> outname, char *outfile, int bin_size, int loci_length, int sample_size, int flag, double up_bound)
{
	char *part;
	char temp_name[255];
	char chrname[10];
	int startsite, count;
	double value;
	FILE * pFile;

	//standard output
	if (flag != 1)
	{
		pFile = fopen(outfile, "w");
		if (pFile != NULL)
		{
			fprintf(pFile, "Chromosome\tStart\tEnd\t");
			for (int i = 0; i < sample_size-1; i++)
			{
				fprintf(pFile, "%s\t", outname[i].c_str());
			}
			fprintf(pFile, "%s\n", outname[sample_size-1].c_str());
			for (int i = 0; i < loci_length; i++)
			{
				fprintf(pFile, "%s\t", select_idx[i]);
				for (int j = 0; j < sample_size-1; j++)
				{
					value = data_out[j][i];
					if (value < 0)
					{
						value = 0;
					}
					else if (value > up_bound)
					{
						value = up_bound;
					}
					fprintf(pFile, "%lf\t", value);   //report original output value, log2(x+1) transformed.
				}

				value = data_out[sample_size-1][i];
				if (value < 0)
				{
					value = 0;
				}
				else if (value > up_bound)
				{
					value = up_bound;
				}
				fprintf(pFile, "%lf\n", value);

			}
			fclose(pFile);
		}
		else
		{
			std::cout << "Error! Cannot write file " << outfile << std::endl;
			return 1;
		}


	}
	else    //WIG output
	{
		std::vector<std::vector<int> > location;
		std::vector<int>::iterator it;
		for (int i = 0; i < 23; i++)
		{
			location.push_back(std::vector<int>());
		}

		for (int i = 0; i < loci_length; i++)
		{
			part = strtok(select_idx[i], "\t\0");
			if (strcmp(part, "chrX") == 0)
			{
				part = strtok(NULL, "\t\0");
				startsite = atoi(part);
				location[22].push_back(startsite);
			}
			else
			{
				for (int j = 0; j < 22; j++)
				{
					sprintf(chrname, "chr%d", j + 1);
					if (strcmp(part, chrname) == 0)
					{
						part = strtok(NULL, "\t\0");
						startsite = atoi(part);
						location[j].push_back(startsite);
						break;
					}
				}
			}
		}

		for (int i = 0; i < sample_size; i++)
		{
			count = 0;
			strcpy(temp_name, outfile);
			strcat(temp_name, ".");
			strcat(temp_name, outname[i].c_str());
			strcat(temp_name, ".wig");
			std::cout << "Writing file " << temp_name << std::endl;
			pFile = fopen(temp_name, "w");
			if (pFile != NULL)
			{
				fprintf(pFile, "track\ttype=wiggle_0\tname=%s\tvisibility=full\tautoScale=off\tmaxHeightPixels=100:50:10\tviewLimits=0.0:100.0\tyLineOnOff=off\n", outname[i].c_str());
				//print chr1-chr22
				for (int j = 0; j < 22; j++)
				{
					fprintf(pFile, "variableStep\tchrom=chr%d\tspan=%d\n", j + 1, bin_size);
					for (it = location[j].begin(); it != location[j].end(); it++)
					{
						value = pow(2, data_out[i][count]) - 1;
						if (value < 0)
						{
							value = 0;
						}
						else if (value > pow(2,up_bound))
						{
							value = pow(2,up_bound);
						}
						count++;
						fprintf(pFile, "%d\t%lf\n", (*it), value);        //report value is limited from 0 to 10000.
					}
				}
				//print chrX
				fprintf(pFile, "variableStep\tchrom=chrX\tspan=%d\n", bin_size);
				for (it = location[22].begin(); it != location[22].end(); it++)
				{
					value = pow(2, data_out[i][count]) - 1;
					if (value < 0)
					{
						value = 0;
					}
					else if (value > pow(2,up_bound))
					{
						value = pow(2,up_bound);
					}
					count++;
					fprintf(pFile, "%d\t%lf\n", (*it), value);
				}

				fclose(pFile);
			}
			else
			{
				std::cout << "Error! Cannot write file " << temp_name << std::endl;
				return 1;
			}

		}
	}
	return 0;
}

// [[Rcpp::export]]
void help_info()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "Standard output: BIRD_predict -b model_file.bin -i input_file.txt -o output_file.txt" << std::endl;
	std::cout << "Standard output will save a matrix contained all predited value in log scale (log2(x+1) transformed)." << std::endl;
	std::cout << "WIG output: BIRD_predict -b model_file.bin -i input_file.txt -o output_name -w" << std::endl;
	std::cout << "WIG output will save each sample as a WIG file." << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-b   Specify library file. If not sepecified,the program will search for model_file.bin in the current directory." << std::endl;
	std::cout << "-i   Specify input file (gene expression obtained from GeneBASE)." << std::endl;
	std::cout << "-o   Specify output file." << std::endl;
	std::cout << "-u   Set upper bound for predicted values (default:14)." << std::endl;
	std::cout << "-w   Output WIG file for each sample." << std::endl;
	std::cout << "-l   Use locus-level model for prediction." << std::endl;
    std::cout << "-e   Use exact id match for matching the gene expression data." << std::endl;
	return;
}

// [[Rcpp::export]]
int main(int argc, char *argv[])
{
	char infile[255];
	char outfile[255];
    char outfile_expr[255];
	char libfile[255]="./model_file.bin";
	int write_flag = 0;
	int locus_model = 0;
	double up_bound = 14;
	int opt;
    int match_mode = 0;
    std::string sep_str = ".";

	if(argc == 1)
	{
		help_info();
		return 1;
	}

	while ((opt = getopt(argc,argv,"b:i:o:u:whle")) != EOF)
	{
		switch(opt)
		{
			case 'b':
				strcpy(libfile, optarg);
				break;
			case 'i':
				strcpy(infile, optarg);
				break;
			case 'o':
				strcpy(outfile, optarg);
				break;
			case 'u':
				up_bound = atof(optarg);
				break;
			case 'w':
				write_flag = 1;
				break;
			case 'h':
				help_info();
				return 1;
			case 'l':
				locus_model = 1;
				break;
            case 'e':
                match_mode = 1;
                break;
			case '?':
				std::cout << "Please input the correct parameters." << std::endl;
				std::cout << "Standard output: BIRD_predict -b model_file.bin -i input_file.txt -o output_file.txt" << std::endl;
				std::cout << "WIG output: BIRD_predict -b model_file.bin -i input_file.txt -o output_name -w" << std::endl;
				return 1;
			default:
				abort();
		}
	}

	//output parameters
	std::cout << "Using model file: " << libfile << std::endl;
	std::cout << "Input file name: " << infile << std::endl;
	std::cout << "Output file name: " << outfile << std::endl;
	std::cout << "Using upper bound: " << up_bound << std::endl;
	if(write_flag)
	{
		std::cout << "Output mode: wig file" << std::endl;
	}
	else
	{
		std::cout << "Output mode: data matrix" << std::endl;
	}

	if(locus_model)
	{
		std::cout << "Prediction mode: locus-level model" << std::endl;
	}
	else
	{
		std::cout << "Prediction mode: full model" << std::endl;
	}

	Exondata exonin;
	Exondata *indata = &exonin;
	int predictor_size, sample_size, cluster_size, var_size, loci_size, bin_size, DH_num1, DH_num2, DH_num3;

	if (ReadPar(libfile, loci_size, predictor_size, cluster_size, bin_size, var_size, DH_num1, DH_num2, DH_num3))
	{
		return 1;
	}

	double *quantile_in = new double[predictor_size];
	int *cluster_idx = new int[predictor_size];
	double *exon_mean = new double[predictor_size];
	double *exon_sd = new double[predictor_size];
	double *DNase_mean = new double[loci_size];
	double *DNase_sd = new double[loci_size];

	double **coef = new double *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		coef[j] = new double[loci_size]();
	}

	int **pre_idx = new int *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		pre_idx[j] = new int[loci_size]();
	}

	char **select_loci = new char *[loci_size];
	for (int j = 0; j < loci_size; j++)
	{
		select_loci[j] = new char[30];
	}

	char **TC_id = new char *[predictor_size];
	for (int j = 0; j < predictor_size; j++)
	{
		TC_id[j] = new char[30];
	}

	//DH cluster model
	double **dis_matrix = new double *[3];
	for (int j = 0; j < 3; j++)
	{
		dis_matrix[j] = new double[loci_size]();
	}

	int **DH_cluster = new int *[3];
	for (int j = 0; j < 3; j++)
	{
		DH_cluster[j] = new int[loci_size]();
	}

	double **DH_coef1 = new double *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_coef1[j] = new double[DH_num1]();
	}

	double **DH_coef2 = new double *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_coef2[j] = new double[DH_num2]();
	}

	double **DH_coef3 = new double *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_coef3[j] = new double[DH_num3]();
	}

	int **DH_pre_idx1 = new int *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_pre_idx1[j] = new int[DH_num1]();
	}

	int **DH_pre_idx2 = new int *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_pre_idx2[j] = new int[DH_num2]();
	}

	int **DH_pre_idx3 = new int *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		DH_pre_idx3[j] = new int[DH_num3]();
	}

	//read in modle file
	if (ReadinModel(libfile, quantile_in, exon_mean, exon_sd, coef, DNase_mean, DNase_sd, pre_idx, TC_id, cluster_idx, select_loci, predictor_size, var_size, loci_size, dis_matrix, DH_cluster, DH_coef1, DH_coef2, DH_coef3, DH_pre_idx1, DH_pre_idx2, DH_pre_idx3, DH_num1, DH_num2, DH_num3))
	{
		return 1;
	}
	//read in gene expression data
	if (ReadinExon(infile, indata))
	{
		return 1;
	}

    int *match_idx = new int[predictor_size];
    
    //match and check gene expression data
    std::cout << "Matching the input gene expression data." << std::endl;
    if(match_mode == 1)
    {
        //use exact match
        if (MatchExon(TC_id, indata, predictor_size, match_idx))
        {
            std::cout << "Gene expression data format incorrect: No gene id matches the library file." << std::endl;
            std::cout << "Please check sample input file for reference." << std::endl;
            return 1;
        }
    }
    else
    {
        //match id before separate character
        if (MatchExon_sep(TC_id, indata, predictor_size, match_idx, sep_str))
        {
            std::cout << "Gene expression data format incorrect: No gene id matches the library file." << std::endl;
            std::cout << "Please check sample input file for reference." << std::endl;
            return 1;
        }
    }
    
    std::cout << "Processing data..." << std::endl;
    
	sample_size = (int)exonin.sample_name.size();

	double ** data_norm = new double *[sample_size];
	double ** data_mean = new double *[sample_size];
	double ** output = new double *[sample_size];
	double ** DH_pre1 = new double *[sample_size];
	double ** DH_pre2 = new double *[sample_size];
	double ** DH_pre3 = new double *[sample_size];

	for (int i = 0; i < sample_size; i++)
	{
		data_norm[i] = new double[predictor_size]();
		data_mean[i] = new double[cluster_size]();
		output[i] = new double[loci_size]();
		DH_pre1[i] = new double[DH_num1]();
		DH_pre2[i] = new double[DH_num2]();
		DH_pre3[i] = new double[DH_num3]();
	}

	for (int i = 0; i < sample_size; i++)
	{
		for (int j = 0; j < predictor_size; j++)
		{
            if(match_idx[j] != -1)
            {
                data_norm[i][j] = log2(exonin.data[match_idx[j]][i] + 1);
            }
            else
            {
                data_norm[i][j] = 0;
            }
		}
	}
    
    strcpy(outfile_expr, infile);
    strcat(outfile_expr, ".match.txt");
    WriteExpr(data_norm, exonin.sample_name, TC_id, outfile_expr, predictor_size, sample_size);
    
    for (int i = 0; i < sample_size; i++)
    {
        QuantileNorm(data_norm[i], quantile_in, predictor_size);
    }

	StandardizeRow(data_norm, exon_mean, exon_sd, predictor_size, sample_size); //standardize gene expression data
	ClusterMean(data_norm, data_mean, cluster_idx, predictor_size, cluster_size, sample_size); //get gene expression cluster mean;

	if(locus_model != 1)
	{
		//Locus level prediction
		Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);

		//DH cluster level prediction
		Regression(data_mean, DH_pre1, DH_coef1, DH_pre_idx1, var_size, DH_num1, sample_size);
		Regression(data_mean, DH_pre2, DH_coef2, DH_pre_idx2, var_size, DH_num2, sample_size);
		Regression(data_mean, DH_pre3, DH_coef3, DH_pre_idx3, var_size, DH_num3, sample_size);

		//Combine result from locus level and DH cluster level
		ModelAverage(output, DH_pre1, DH_pre2, DH_pre3, dis_matrix, DH_cluster, loci_size, sample_size);
		std::cout << "Using full model for prediction." << std::endl;
	}
	else
	{
		//Locus level prediction
		Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);
		std::cout << "Using locus-level model for prediction." << std::endl;
	}

  	//convert predicted value back to original scale
	StandardizeRow_r(output, DNase_mean, DNase_sd, loci_size, sample_size);

  	//write output file
	std::cout << "Writing output file..." << std::endl;
	if (WriteWIG(output, select_loci, exonin.sample_name, outfile, bin_size, loci_size, sample_size, write_flag, up_bound))
	{
		return 1;
	}

	//release memory.
	ReleaseExondata(exonin);
	delete[] quantile_in;
	delete[] cluster_idx;
	for (int i = 0; i < sample_size; i++)
	{
		delete[] data_norm[i];
		delete[] data_mean[i];
		delete[] output[i];
		delete[] DH_pre1[i];
		delete[] DH_pre2[i];
		delete[] DH_pre3[i];
	}
	for (int i = 0; i < var_size; i++)
	{
		delete[] coef[i];
		delete[] pre_idx[i];
		delete[] DH_coef1[i];
		delete[] DH_coef2[i];
		delete[] DH_coef3[i];
		delete[] DH_pre_idx1[i];
		delete[] DH_pre_idx2[i];
		delete[] DH_pre_idx3[i];
	}
	for (int i = 0; i < 3; i++)
	{
		delete[] dis_matrix[i];
		delete[] DH_cluster[i];
	}
	for (int i = 0; i < predictor_size; i++)
	{
		delete[] TC_id[i];
	}
	for (int i = 0; i < loci_size; i++)
	{
		delete[] select_loci[i];
	}
	delete[] data_norm;
	delete[] data_mean;
	delete[] output;
	delete[] DH_pre1;
	delete[] DH_pre2;
	delete[] DH_pre3;
	delete[] coef;
	delete[] pre_idx;
	delete[] DH_coef1;
	delete[] DH_coef2;
	delete[] DH_coef3;
	delete[] DH_pre_idx1;
	delete[] DH_pre_idx2;
	delete[] DH_pre_idx3;
	delete[] TC_id;
	delete[] select_loci;
	delete[] dis_matrix;
	delete[] DH_cluster;
	return 0;
}
