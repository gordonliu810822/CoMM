//
//  plinkfun.cpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//
#include <RcppArmadillo.h>
#include "plinkfun.hpp"
#include <stdio.h>
#include <math.h>
#include <bitset>
#include <boost/algorithm/string.hpp>

#define MAX_LEN 20
int getLineNum(string filename){
    FILE *pf = fopen(filename.c_str(), "r"); // 打开文件
    char buf[10000];
    int lineCnt = 0;
    if (!pf) // 判断是否打开成功
        return -1;
    while (fgets(buf, 10000, pf)) // fgets循环读取，直到文件最后，才会返回NULL
        lineCnt++; // 累计行数
    fclose(pf);
    return lineCnt;
}

void getFourGentype2(char* geno, std::bitset<8> bits){
    int idx = 0;
    for (int j=0; j < 8; j = j + 2) {
        if(bits[j] && bits[j+1]){
            geno[idx] = 0;
        }else if(!bits[j] && !bits[j+1]){
            geno[idx] = 2;
        }else if(!bits[j] && bits[j+1]){
            geno[idx] = 1;
        }else if(bits[j] && !bits[j+1]){
            geno[idx] = 3;
        }
        idx++;
    }
}

void getFourGentype(unsigned* geno, std::bitset<8> bits){
    int idx = 0;
    for (int j=0; j < 8; j = j + 2) {
        if(bits[j] && bits[j+1]){
            geno[idx] = 0;
        }else if(!bits[j] && !bits[j+1]){
            geno[idx] = 2;
        }else if(!bits[j] && bits[j+1]){
            geno[idx] = 1;
        }else if(bits[j] && !bits[j+1]){
            geno[idx] = 3;
        }
        idx++;
    }
}

void readPlink(string stringname, long int N, long int P, unsigned* X){

   // string stringname = dir + dataname;
    FILE *fp;
    unsigned char buff[3];
    string bedfile = stringname +".bed";
    fp = fopen(bedfile.c_str(), "rb");
    if (!fp) return;
    fread(buff, sizeof(char), 3, fp);

    std::bitset<8> magic1(buff[0]);
    std::bitset<8> magic2(buff[1]);
    std::bitset<8> mode0(buff[2]);

    if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
     //   cout <<"Error Identifier of plink binary file" << endl;
    }

    unsigned long mode =  mode0.to_ulong();
    if(mode == 0){
        printf ("individual-Major Order:improper type of plink file");
        exit (EXIT_FAILURE);
    }
    //     cout << "SNP-Major Order" << endl;
    // }else if(mode == 0){
    //    cout << "individual-Major Order" << endl;
    // }
    // X = new int[N*P];
    //cout << N << "x" << P << endl;
    //cout << sizeof(*X) << endl;
    long n = 0;
    long long charNum = ceil(N*1.0/4)*10000;
    long long leftGenoNum = ceil(N*1.0/4)*P;
    long nblock = ceil(N*1.0/4);
    long nSNP = 0;
    //cout << "nblock: " << nblock << endl;
    //cout << "leftGenoNum: " << leftGenoNum << endl;
    while (!feof(fp)) {
        if(leftGenoNum <= 0)
            break;
        if(leftGenoNum <= charNum){
            charNum  = leftGenoNum;
        }
        char* genotype = new char[charNum];
        fread(genotype, sizeof(char), charNum, fp);
        unsigned* geno = new unsigned[4];
        long nSNPc = long(charNum / nblock); //number of SNPs of this iteration

        //cout << n << "-th line: ";
        //cout << "nSNPc: "<< nSNPc << endl;
        // cout << sizeof(int) << endl;
        long long idx = 0;
        for (long i=0; i < nSNPc; i++) {
			
            /*if ( n >= 4 ){
                cout << i << "-th snp" << endl;
				if (i == 0){
					//cout << "break 1 ... " << nSNP << ";" << N << ";" << leftGenoNum << ";" << idx << endl;
				}
            }*/
            for(long j=0; j < nblock - 1; j++){
				long long indx = (long long)(nSNP) * (long long)(N) + (long long)(j*4);

                // if ( n == 3 && i == 9999  && j == nblock - 2){
				/*if ( n == 3 && j == nblock - 2 && i > 5000 ){
					// long long indxt = nSNP * N + j*4;
					// cout << "break 2 ... " << endl;
					cout << "break 1 ... "<< n << "-th line, " << i << "-th snp: " << N <<";" << j  << ";" << 
						indx << ";" << leftGenoNum << ";" << idx << endl;
                }*/

				std::bitset<8> bits(genotype[idx]);
				getFourGentype(geno,bits);
				memcpy(X + indx, geno, 4*sizeof(unsigned));
			

                idx++;
                leftGenoNum -= 1;
            }
			//if ( n == 4  && i == 0 ){
			//	cout << "break 4 ... " << endl;
			//}
			// cout << "break 1 ... " << endl;
            long left = N - (nblock - 1)*4;
            std::bitset<8> bits(genotype[idx]);
            getFourGentype(geno,bits);

			long long indx2 = (long long)(nSNP) * (long long)(N) + (long long)(nblock - 1)*4;
			long long indx3 = left*sizeof(unsigned);
			/*if ( n == 3 && i > 5000 ){
				cout << "break 2 ... " << n << "-th line, " << i << "-th snp: " << nSNP << ";" << 
					N << ";" << nblock << ";" << indx2 << ";" << indx3 << endl;
			}*/
            memcpy(X + indx2, geno, indx3);
            idx++;
            leftGenoNum -= 1;
            nSNP ++;
        }

        delete[] geno;
        delete[] genotype;
        n++;
    //    cout <<n << " processing"<<endl;
    }


}

void readPlink2(string stringname,int N, int P, char* X){

    FILE *fp;
    unsigned char buff[3];
    string bedfile = stringname +".bed";
    fp = fopen(bedfile.c_str(), "rb");
    if (!fp) return;
    fread(buff, sizeof(char), 3, fp);

    std::bitset<8> magic1(buff[0]);
    std::bitset<8> magic2(buff[1]);
    std::bitset<8> mode0(buff[2]);

    if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
     //   cout <<"Error Identifier of plink binary file" << endl;
    }
    unsigned long mode =  mode0.to_ulong();
    if(mode == 0){
        printf ("individual-Major Order:improper type of plink file");
        exit (EXIT_FAILURE);
    }
    long n = 0;
    long long charNum = ceil(N*1.0/4)*10000;
    long long leftGenoNum = ceil(N*1.0/4)*P;
    long nblock = ceil(N*1.0/4);
    long nSNP = 0;
    while (!feof(fp)) {
        if(leftGenoNum <= 0)
            break;
        if(leftGenoNum <= charNum){
            charNum  = leftGenoNum;
        }
        char* genotype = new char[charNum];
        fread(genotype, sizeof(char), charNum, fp);
        char* geno = new char[4];
        long nSNPc = long(charNum / nblock); //number of SNPs of this iteration
        long long idx = 0;
        for (long i=0; i < nSNPc; i++) {
          for(long j=0; j < nblock - 1; j++){
            long long indx = (long long)(nSNP) * (long long)(N) + (long long)(j*4);
            std::bitset<8> bits(genotype[idx]);
            getFourGentype2(geno,bits);
            memcpy(X + indx, geno, 4*sizeof(char));
            idx++;
            leftGenoNum -= 1;
          }
          long left = N - (nblock - 1)*4;
          std::bitset<8> bits(genotype[idx]);
          getFourGentype2(geno,bits);

          long long indx2 = (long long)(nSNP) * (long long)(N) + (long long)(nblock - 1)*4;
          long long indx3 = left*sizeof(char);
          memcpy(X + indx2, geno, indx3);
          idx++;
          leftGenoNum -= 1;
          nSNP ++;
        }

        delete[] geno;
        delete[] genotype;
        n++;
    }
}

mat getSubMat(char* X, int N , int P, uvec indics, double* sub_matrix_double){
    int idx = 0;
    int elem = 0;
    char* start_p = NULL;
    for(int i = 0; i < indics.n_elem; i++){
        start_p = X + indics[i] * N;
        for(int j = 0; j < N; j++){
            elem = start_p[j];
            sub_matrix_double[idx++] = elem;
        }
    }
    int subP = indics.n_elem;
    Mat<double> A(sub_matrix_double, N, subP, false);
    return A;
}


void ReadPlinkBimFile(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
                      IntegerVector chr, IntegerVector bp, NumericVector morgan, int P){
    
    
    FILE *stream;
    /*CharacterVector rsname(P);
     IntegerVector A1(P), A2(P);*/
    
    int ch, b;
    double mor;
    char s[MAX_LEN + 1], efa, nefa;
    
    stream = fopen(stringname.c_str(), "r");
    clock_t t1 = clock();
    
    //int i = 0;
    /* Put in various data. */
    for ( int i = 0; i < P; i++){
        if (i % 500000 == 0 && i != 0){
            cout << i << "-th SNP" << ",";
            cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
        }
        
        fscanf(stream, "%i %s %lf %i %c %c", &ch, &s[0], &mor, &b, &efa, &nefa);
        
        chr(i) = ch;
        rsname(i) = s;
        morgan(i) = mor;
        bp(i) = b;
        A1(i) = (int)efa;
        A2(i) = (int)nefa;
        
    }
    
    
}

CharacterVector charv_subset(CharacterVector x, uvec idx){
    CharacterVector v(idx.n_elem);
    for (unsigned int i = 0; i < idx.n_elem; i++){
        v(i) = x(idx(i));
    }
    return v;
}

void ReadPlinkFamFile(std::string stringname, CharacterVector FID, CharacterVector IID, IntegerVector sex,
                      NumericVector pheno, int N)
{
    FILE *stream;
    int gender;
    double phn;
    char fid[MAX_LEN + 1], iid[MAX_LEN + 1], tmp1[MAX_LEN + 1], tmp2[MAX_LEN + 1];
    
    stream = fopen(stringname.c_str(), "r");
    
    for (int i = 0; i < N; i++){
        fscanf(stream, "%s %s %s %s %i %lf", &fid, &iid, &tmp1, &tmp2, &gender, &phn);
        
        FID(i) = fid;
        IID(i) = iid;
        sex(i) = gender;
        pheno(i) = phn;
        
    }
}

void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID,
                       NumericVector pheno, int nrows, int whCol){
    
    std::ifstream myfile(stringname.c_str());
    std::string line;
    
    clock_t t1 = clock();
    
    int nrow_ind = 0;
    vector <string> tmp;
    
    if (myfile.is_open()){
        while (nrow_ind < nrows){
            if (nrow_ind % 5000 == 0 && nrow_ind != 0){
                cout << nrow_ind << "-th individual " << ",";
                cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
            }
            
            getline(myfile, line);
            boost::split(tmp, line, boost::is_any_of(" \t *"));
            
            FID(nrow_ind) = tmp[0];
            IID(nrow_ind) = tmp[1];
            pheno(nrow_ind) = atof(tmp[4 + whCol].c_str());
            
            nrow_ind++;
        }
    }
}
