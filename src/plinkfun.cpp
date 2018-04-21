//
//  plinkfun.cpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//
#include <RcppArmadillo.h>
//#include <Rcpp.h>
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
/*void readPlink(string stringname,int N, int P, unsigned* X){

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
			
            //if ( n >= 4 ){
            //    cout << i << "-th snp" << endl;
			//	if (i == 0){
			//		//cout << "break 1 ... " << nSNP << ";" << N << ";" << leftGenoNum << ";" << idx << endl;
			//	}
            //}
            for(long j=0; j < nblock - 1; j++){
				long long indx = (long long)(nSNP) * (long long)(N) + (long long)(j*4);

                // if ( n == 3 && i == 9999  && j == nblock - 2){
				//if ( n == 3 && j == nblock - 2 && i > 5000 ){
					// long long indxt = nSNP * N + j*4;
					// cout << "break 2 ... " << endl;
				//	cout << "break 1 ... "<< n << "-th line, " << i << "-th snp: " << N <<";" << j  << ";" << 
				//		indx << ";" << leftGenoNum << ";" << idx << endl;
                //}

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
			//if ( n == 3 && i > 5000 ){
			//	cout << "break 2 ... " << n << "-th line, " << i << "-th snp: " << nSNP << ";" << 
			//		N << ";" << nblock << ";" << indx2 << ";" << indx3 << endl;
			//}
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


}*/


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
			// sex(nrow_ind) = atoi(tmp[4].c_str());
			pheno(nrow_ind) = atof(tmp[4 + whCol].c_str());

			//cout << "value: " << tmp[0] << ";" << tmp[1] << ";" << tmp[2] << ";" << tmp[3] << ";" << tmp[4]
			//	<< ";" << tmp[5] << ";" << tmp[6] << endl;
			nrow_ind++;
		}
	}
}

List match_SNPs(std::string stringname1, std::string stringname2){

	// stringname1: prefix for plink file 1; start working on plink file 1 (1000 G with expression data)
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	//cout << "break 1: " << bimfile1  << endl;
	int P1 = getLineNum(bimfile1);
	cout << "## Number of SNPs (plink file 1):" << P1 << endl;

	IntegerVector A1_1(P1), A2_1(P1);
	CharacterVector rsname_1(P1);
	IntegerVector chr_1(P1), bp_1(P1);
	NumericVector morgan_1(P1);

	ReadPlinkBimFile(bimfile1, A1_1, A2_1, rsname_1, chr_1, bp_1, morgan_1, P1);
	// cout << rsname_1(0) << ";" << A1_1(0) << ";" << A2_1(0) << endl;

	// stringname2: prefix for plink file 2; start working on plink file 2 (GWAS data with trait)
	string bimfile2 = stringname2;
	bimfile2 += ".bim";

	int P2 = getLineNum(bimfile2);
	cout << "## Number of SNPs (plink file 2):" << P2 << endl;

	IntegerVector A1_2(P2), A2_2(P2);
	CharacterVector rsname_2(P2);
	IntegerVector chr_2(P2), bp_2(P2);
	NumericVector morgan_2(P2);

	ReadPlinkBimFile(bimfile2, A1_2, A2_2, rsname_2, chr_2, bp_2, morgan_2, P2);
	// cout << rsname_2(0) << ";" << chr_2(0) << ";" << bp_2(0) << ";" << A1_2(0) << ";" << A2_2(0) << endl;

	// mathcing panel SNPs in file 1 and file 2 with correction for direction of ref allele
	// rsname in both files in the order of the first file.
	CharacterVector rs_inter = intersect(rsname_1, rsname_2);
	IntegerVector idxin = match(rsname_1, rs_inter); //index for SNPs in file 1
	CharacterVector rsname_4use = rsname_1[Rcpp::is_na(idxin) == false];
	IntegerVector chr_4use = chr_1[Rcpp::is_na(idxin) == false];
	IntegerVector bp_4use = bp_1[Rcpp::is_na(idxin) == false];

	// match snps (rsname_4use; rsname_2: file 2)
	IntegerVector idxin2 = match(rsname_4use, rsname_2);  //index for SNPs in file 2
	// match snps (rsname_4use; rsname_1: file 1)
	IntegerVector idxin1 = match(rsname_4use, rsname_1);  //index for SNPs in file 1

	/*IntegerVector chr_4use_tmp = chr_2[idxin2];
	IntegerVector bp_4use_tmp = bp_2[idxin2];

	vec idxtmp = as<vec>(chr_4use) - as<vec>(chr_4use_tmp);
	cout << "check the quality: " << sum(idxtmp) << endl;

	cout << "Size of matched SNPs: " << rsname_4use.size() << endl;*/

	// convert ascii letter to uvec and work on overlapped SNPs
	uvec idx = as<uvec>(idxin1) -1;
	uvec tmp1 = as<uvec>(A1_1);
	uvec A1_1_ = tmp1.elem(idx);
	tmp1 = as<uvec>(A2_1);
	uvec A2_1_ = tmp1.elem(idx);

	idx = as<uvec>(idxin2) -1;
	uvec tmp2 = as<uvec>(A1_2);
	uvec A1_2_ = tmp2.elem(idx);
	tmp2 = as<uvec>(A2_2);
	uvec A2_2_ = tmp2.elem(idx);

	//ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
	/*//replace lower letter to upper letter in A1_1_ and A2_1_;
	idx = find(A1_1_ > 85);
	A1_1_.elem(idx) = A1_1_.elem(idx) - 32;
	idx = find(A2_1_ > 85);
	A2_1_.elem(idx) = A2_1_.elem(idx) - 32;

	idx = find(A1_2_ > 85);
	A1_2_.elem(idx) = A1_2_.elem(idx) - 32;
	idx = find(A2_2_ > 85);
	A2_2_.elem(idx) = A2_2_.elem(idx) - 32;*/

	//compare A1_1_ A1_2_, A2_1_ A2_2_
	//A1_1_: replace T with A,replace G with C
	idx = find(A1_1_ == 84);
	uvec idxrepl1(idx.n_elem);
	idxrepl1.fill(65);
	A1_1_.elem(idx) = idxrepl1;

	idx = find(A1_1_ == 71);
	uvec idxrepl2(idx.n_elem);
	idxrepl2.fill(67);
	A1_1_.elem(idx) = idxrepl2;

	//A1_2_: replace T with A,replace G with C
	idx = find(A1_2_ == 84);
	uvec idxrepl3(idx.n_elem);
	idxrepl3.fill(65);
	A1_2_.elem(idx) = idxrepl3;

	idx = find(A1_2_ == 71);
	uvec idxrepl4(idx.n_elem);
	idxrepl4.fill(67);
	A1_2_.elem(idx) = idxrepl4;

	//A2_1_: replace T with A,replace G with C
	idx = find(A2_1_ == 84);
	uvec idxrepl5(idx.n_elem);
	idxrepl5.fill(65);
	A2_1_.elem(idx) = idxrepl5;

	idx = find(A2_1_ == 71);
	uvec idxrepl6(idx.n_elem);
	idxrepl6.fill(67);
	A2_1_.elem(idx) = idxrepl6;

	//A2_2_: replace T with A,replace G with C
	idx = find(A2_2_ == 84);
	uvec idxrepl7(idx.n_elem);
	idxrepl7.fill(65);
	A2_2_.elem(idx) = idxrepl7;

	idx = find(A2_2_ == 71);
	uvec idxrepl8(idx.n_elem);
	idxrepl8.fill(67);
	A2_2_.elem(idx) = idxrepl8;

	// remove index;
	idx = find((A1_1_ + A2_1_) == (A1_2_ + A2_2_));
	uvec tmp3 = as<uvec>(idxin2);
	uvec idxin22 = tmp3.elem(idx);
	tmp3 = as<uvec>(idxin1);
	uvec idxin11 = tmp3.elem(idx);

	uvec A1_1_r = A1_1_.elem(idx), A2_1_r = A2_1_.elem(idx);
	uvec A1_2_r = A1_2_.elem(idx), A2_2_r = A2_2_.elem(idx);

	uvec chr_tmp = as<uvec>(chr_4use);
	uvec chr_4use_r = chr_tmp.elem(idx);
	uvec bp_tmp = as<uvec>(bp_4use);
	uvec bp_4use_r = bp_tmp.elem(idx);

	CharacterVector rsname_4use_r = charv_subset(rsname_4use, idx);

	vec ind(A1_1_r.n_elem);
	idx = find(A1_1_r == A1_2_r);
	ind.ones();
	ind = -ind;
	ind.elem(idx).ones();

	cout << "Number of matched SNPs (plink file 1 and 2):" << ind.size() << endl;
	//cout << "direction (compare with trait_1000Q_match.R): " << sum(ind) << endl; //compare sum(ind) in R: Height_1000Q_match.R
	//cout << "Size of matched SNPs (remove ambiguous SNPs): " << A1_1_r.n_elem << endl;

	List output = List::create(Rcpp::Named("chr_4use_r") = chr_4use_r,
		Rcpp::Named("A1_1_r") = A1_1_r,
		Rcpp::Named("A1_2_r") = A1_2_r,
		Rcpp::Named("A2_1_r") = A2_1_r,
		Rcpp::Named("A2_2_r") = A2_2_r,
		Rcpp::Named("bp_4use_r") = bp_4use_r,
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("indicator") = ind,
		Rcpp::Named("idxinFile1") = idxin11,
		Rcpp::Named("idxinFile2") = idxin22);

	return output;
}



