rm(list = ls())
library(CoMM)
ls("package:CoMM")

chr = 22

bw = 500000
k=5

  
#file1 = paste("D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\Debug\\1000G.EUR.QC.",chr,sep="");
file1 = paste("/Users/nus/GoogleDrive/Work/Multiple_Trait/TWAS/Debug/1000G.EUR.QC.",chr,sep="");
#file1 = "/Volumes/New Volume/Data/GWAS/1000genomeEuro/1000G_EUR_Phase3_plink/1000G.EUR.QC.ph3.geno001.mind005.hwe0001maf001";
#file1 = paste("M:\\Data\\GWAS\\1000genomeEuro\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.hm3.ph3.MAF001HWE00001",sep="");
#file2 = "D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\Debug\\NFBC_filter_mph10.22";
#file2 = "/Volumes/New Volume/Data/GWAS/NFBC/NFBC_filter_mph10";
#file2 = "M:\\Data\\GWAS\\NFBC\\NFBC_filter_mph10";
file2 = "/Users/nus/GoogleDrive/Work/Multiple_Trait/TWAS/Debug/NFBC_filter_mph10.22";
#file3 = "D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\Debug\\Geuvadis_gene_expression_qn.txt";
file3 = "/Users/nus/GoogleDrive/Work/Multiple_Trait/TWAS/Debug/Geuvadis_gene_expression_qn.txt";
file4 = "";#paste("G:\\Data\\GWAS\\1000genomeEuro\\1000G_EUR_Phase3_plink\\pc4_1000G.EUR.QC.",chr,".txt",sep="");
file5 = "";#"G:\\Data\\GWAS\\NFBC\\pc4_NFBC_filter_mph10.txt";

whichPheno = k;

start_time <- Sys.time() 
fmm = CoMM_testing_run_mt(file1,file2,file3, file4,file5, whichPheno, bw,6);
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time() 
fm = CoMM_testing_run(file1,file2,file3, file4,file5, whichPheno, bw);
end_time <- Sys.time()
end_time - start_time

diff=fm$out_param0-t(fmm$out_param0);
max(abs(diff))
sum(diff^2)

sum(fm$out_param0-fmm$out_param)

out=fm$out_param
outt=out[out[,4]!=-99,]
#out=fm$out_param
#exprinfo = fm$expr_info[out[,4]!=-99,]
#sum(fm$genetype1[out[,4]!=-99] == fm$gene_type1)
#sum(fm$genetype2[out[,4]!=-99] == fm$gene_type2)
out=fm$out_param
outm=fmm$out_param
index=1:dim(out)[1]
idx=index[out[,4]!=-99]


## End loading files ... 
20-th gene: sigma2y: 0.99291; sigma2beta: 0.0352435var of y: 1.03092; dim of y: 344-by-1; var of z: 0.999939; dim of z: 5123-by-1; var of X1tmp[,1]: 0.00456621; dim of X1tmp: 344-by-219; var of X2tmp[,1]: 0.00456621; dim of X2tmp: 5123-by-219; beta0:    0.0076
; dim of beta0: 1-by-1; sum of w1tmp:   3.4400e+002
; dim of w2: 5123-by-1; sum of w2:   5.1230e+003









