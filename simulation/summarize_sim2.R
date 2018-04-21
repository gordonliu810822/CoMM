#rm(list = ls())
#setwd("D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\AUDI\\Simulation\\Archive")
#
#L = 1; M = 200;
#n1 = 350; n2 = 5000;
#q1 = 1; q2 = 5;
#nrep = 500;
#
##rej_prop = array(0,dim=c(3,4,4,5,3));
##alpha_est = array(0,dim=c(3,4,4,5,2));
#rhoall = c(0.2,0.5,0.8)#i
#h2yall = c(0.01,0.05,0.1,0.2)#j
#beta_propall = c(0.02, 0.1, 0.2, 1)#k
#h2all = c(0.0, 0.01, 0.02, 0.05)#l
#pvalall = NULL;#matrix(nrow=500*3*4*4*4*4,ncol=5)
#alphaall = NULL;#matrix(nrow=500*3*4*4*4*3,ncol=5)
#rej_prop = NULL;
#for ( i in 1:3){
#    for (j in 1:3){
#        for (k in 1:4){
#            for ( l in 1:4){
#                rho = rhoall[i];
#                h2y = h2yall[j];
#                beta_prop = beta_propall[k];
#                h2 = h2all[l];
#
#                outfile1 <- paste("sim_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
#                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
#                outfile2 <- paste("sim_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
#                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
#
#                pval <- matrix(unlist(read.table(outfile1,header=F)),ncol=4);
#                alpha <- matrix(unlist(read.table(outfile2,header=F)),ncol=3);
#                #rej_prop[i,j,k,l,] = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T)
#                #alpha_est[i,j,k,l,] = apply(alpha, 2 , mean,na.rm=T)
#                pvalall = rbind(pvalall,data.frame( as.numeric(as.character(pval)),rep(rho,500*4),
#                        rep(h2y,500*4),rep(beta_prop,500*4), rep(h2,500*4),
#                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"),each=500),
#                        rep(M,500*4), rep(n2,500*4)))
#                alphaall= rbind(alphaall,data.frame( as.numeric(as.character(alpha)), rep(rho,500*3),
#                        rep(h2y,500*3),rep(beta_prop,500*3), rep(h2,500*3),
#                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet"),each=500),
#                        rep(M,500*3), rep(n2,500*3)))
#                tmp_rej = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T);
#                rej_prop = rbind(rej_prop,data.frame(cbind(tmp_rej,rep(rho,4),rep(h2y,4),rep(beta_prop,4),
#                        rep(h2,4),c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"), rep(M,4),rep(n2,4))))
#            }
#
#        }
#
#    }
#
#}
#
#pvalall = data.frame(pvalall)
#colnames(pvalall) = c("Pval","RhoX","h2y","beta_prop","h2","Method","M","n2")
#alphaall = data.frame(alphaall)
#colnames(alphaall) = c("Alpha","RhoX","h2y","beta_prop","h2","Method","M","n2")
#
#rej_prop = data.frame(rej_prop);
#colnames(rej_prop) = c("rej_prop","RhoX","h2y","beta_prop","h2","Method","M","n2")
#
#save(pvalall, alphaall, rej_prop, file = "sim1.RData")
#
############################################
##rhoall = c(0.2,0.5,0.8)#i
##h2yall = c(0.01,0.05,0.1,0.2)#j
##beta_propall = c(0.02, 0.1, 0.2, 1)#k
##h2all = c(0.0, 0.01, 0.02, 0.05)#l
#library(ggplot2)
##t1e_pval = pvalall[,c(1,2,3,4,6)]
#
#t1e_pval = pvalall[pvalall$RhoX==0.2,]
#t1e_pval$beta_prop = factor(t1e_pval$beta_prop)
#t1e_pval$h2y = factor(t1e_pval$h2y)
#t1e_pval$h2 = factor(t1e_pval$h2)
#t1e_pval$RhoX = factor(t1e_pval$RhoX)
#t1e_pval$Method <- ordered(t1e_pval$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))
#
#P_pval   <- ggplot(t1e_pval, aes(x = beta_prop, y = -log10(Pval),fill = Method))+ geom_boxplot()+
#geom_hline(yintercept=0.3,colour="orange",linetype="dashed") +
#facet_grid(h2~h2y) +
#theme(legend.position="bottom")
#P_pval
#outfile = gsub("( )", "", paste("h2_005.png"))
#ggsave(outfile, width = 320, height = 240, units = "mm", dpi=600)
#
#
##P_pval   <- ggplot(t1e_pval,aes(x=n2,y=-log10(Pval),color=Method)) + geom_boxplot()  + geom_hline(yintercept=1,colour="orange",linetype="dashed") + facet_grid(snrz~snry,labeller = label_both,scales = "free_y") + theme(legend.position="bottom")
##P_pval
#
##t1e_rej = rej_prop[rej_prop$h2==0.05,c(1,2,3,4,6)]
##t1e_rej$beta_prop = factor(t1e_rej$beta_prop)
##t1e_rej$h2y = factor(t1e_rej$h2y)
##t1e_rej$RhoX = factor(t1e_rej$RhoX)
##t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))
##
##ggplot(t1e_rej, aes(x = h2y, y = rej_prop,color = Method))+
##geom_bar(stat="identity",fill="white", position=position_dodge(preserve = "single"))+
##facet_grid(RhoX~beta_prop) +
##theme(legend.position="bottom")
#for (rho in c(0.2,0.5,0.8)){
#for (beta_prop in c(0.02, 0.1, 0.2, 1)){
#for (n2 in 5000){
#t1e_rej = rej_prop[rej_prop$RhoX==rho&rej_prop$beta_prop==beta_prop&rej_prop$n2==n2,]
##t1e_rej = rej_prop[rej_prop$RhoX==0.2&rej_prop$beta_prop==1&t1e_rej$h2==0&t1e_rej$h2y==0.4&t1e_rej$n2==5000,]
#t1e_rej$h2 = factor(t1e_rej$h2)
#t1e_rej$h2y = factor(t1e_rej$h2y)
#t1e_rej$n2 = factor(t1e_rej$n2)
#t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))
#t1e_rej$rej_prop = as.numeric(as.character((t1e_rej$rej_prop)))
#
##rej_propu = filter(rej_prop, rank(rej_prop) > 1)
#ggplot(t1e_rej, aes(x = Method, y = rej_prop,fill = Method))+
#geom_bar(stat="identity")+
#facet_grid(h2y~h2) + scale_y_continuous(breaks=seq(0,1,0.05))+ 
#theme(panel.grid.major = element_blank(),legend.position="bottom",axis.ticks.y = element_blank())
#
#outfile = paste("barplot_rho_",gsub('[.]','',as.character(rho)),
#"_bp_",gsub('[.]','',as.character(beta_prop)),"_n2_",n2,".png",sep="")
#ggsave(outfile, width = 320, height = 240, units = "mm", dpi=600)
#}}}
#


## 0907
rm(list = ls())
setwd("D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\AUDI\\Simulation\\sim")

L = 1; M = 100;
n1 = 250; #n2 = 5000;
q1 = 1; q2 = 1;

#rej_prop = array(0,dim=c(3,4,4,5,3));
#alpha_est = array(0,dim=c(3,4,4,5,2));
nrep = 500;
rhoall = c(-0.8,-0.5,-0.2,0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0, 0.001, 0.002, 0.003);
n2all = c(4000,8000)

pvalall = NULL;#matrix(nrow=500*3*4*4*4*4,ncol=5)
alphaall = NULL;#matrix(nrow=500*3*4*4*4*3,ncol=5)
rej_prop = NULL;
na_prop = NULL;
for ( i in 1:6){
    for (j in 1:5){
        for (k in 1:6){
            for ( l in 1:4){
              for (nn in 1:2){
                rho = rhoall[i];
                h2y = h2yall[j];
                beta_prop = beta_propall[k];
                h2 = h2all[l];
                n2 = n2all[nn]
		
		if (h2 == 0) {
                    outfile1 <- paste("sim_t1e_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                    outfile2 <- paste("sim_t1e_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
		}
                if (h2 != 0) {
                    outfile1 <- paste("sim_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                    outfile2 <- paste("sim_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
	        }
                pval <- matrix(unlist(read.table(outfile1,header=F)),ncol=5);
                alpha <- matrix(unlist(read.table(outfile2,header=F)),ncol=4);
                pval[pval==0]=1e-50
                pval[is.na(pval)] = 1;
     
                pvalall = rbind(pvalall,data.frame( as.numeric(as.character(pval)),rep(rho,nrep*5),
                        rep(h2y,nrep*5),rep(beta_prop,nrep*5), rep(h2,nrep*5),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT","2-stage:AUDI"),each=nrep),
                        rep(M,nrep*5), rep(n2,nrep*5)))
                alphaall= rbind(alphaall,data.frame( as.numeric(as.character(alpha)), rep(rho,nrep*4),
                        rep(h2y,nrep*4),rep(beta_prop,nrep*4), rep(h2,nrep*4),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","2-stage:AUDI"),each=nrep),
                        rep(M,nrep*4), rep(n2,nrep*4)))
                tmp_rej = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T);
                sd_rej =  sqrt(tmp_rej*(1-tmp_rej)/nrep)
                tmp_na_rej = apply(is.na(pval) ,2,sum) / nrep;
                rej_prop = rbind(rej_prop,data.frame(cbind(tmp_rej,tmp_na_rej,sd_rej,rep(rho,5),rep(h2y,5),
                        rep(beta_prop,5),
                        rep(h2,5),
                        c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT","2-stage:AUDI"), rep(M,5),rep(n2,5))))
                #na_prop = rbind(na_prop, data.frame(cbind(tmp_na_rej,rep(rho,4),rep(h2y,4),rep(beta_prop,4),
                #        rep(h2,4),c("AUDI","2-stage:Ridge","2-stage:Lasso","SKAT"), rep(M,4),rep(n2,4))))
            }
            }

        }

    }

}

pvalall = data.frame(pvalall)
colnames(pvalall) = c("Pval","RhoX","h2y","beta_prop","h2","Method","M","n2")
alphaall = data.frame(alphaall)
colnames(alphaall) = c("Alpha","RhoX","h2y","beta_prop","h2","Method","M","n2")
rej_prop = data.frame(rej_prop);
colnames(rej_prop) = c("rej_prop","na_prop","sd_rej","RhoX","h2y","beta_prop","h2","Method","M","n2")
#na_prop = data.frame(na_prop);
#colnames(na_prop) = c("na_prop","RhoX","h2y","beta_prop","h2","Method","M","n2")

dat = list();
dat[[1]] = pvalall;
dat[[2]] = alphaall;
dat[[3]] = rej_prop;

save(dat, file = "sim3_n1_250.RData")

##boxplot
#library(ggplot2)
#
#for (rho in c(0.2,0.5,0.8)){
#for (n2 in c(10000)){
#t1e_pval = pvalall[pvalall$RhoX==rho&pvalall$n2==n2,]
#t1e_pval$beta_prop = factor(t1e_pval$beta_prop)
#t1e_pval$h2 = factor(t1e_pval$h2)
#t1e_pval$h2y = factor(t1e_pval$h2y)
#t1e_pval$Method <- ordered(t1e_pval$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Lasso","SKAT"))
#t1e_pval$n2 = factor(t1e_pval$n2);
#
#
#P_pval   <- ggplot(t1e_pval, aes(x = beta_prop, y = -log10(Pval),fill = Method))+ geom_boxplot()+
#geom_hline(yintercept=0.3,colour="orange",linetype="dashed") +
#facet_grid(h2~h2y,labeller = label_both,scales = "free_y") +
#theme(legend.position="bottom",axis.text=element_text(size=18),
#axis.title.x=element_text(size=32),axis.title.y=element_text(size=32),
#legend.title=element_text(size=20), legend.text=element_text(size=20))
#P_pval
#outfile = paste("boxplot_rho_",gsub('[.]','',as.character(rho)),
#"_n2_",gsub('[.]','',as.character(n2)),".png",sep="")
#ggsave(outfile, width = 320, height = 240, units = "mm", dpi=600)
#}}

############################################
##barplot
#dat_rej1 = rej_prop[,2:9]
#colnames(dat_rej1) = c("Percent","RhoX","h2y","beta_prop","h2","Method","M","n2")
#dat_rej1$PropC = rep("NA", 648);
#dat_rej2 = rej_prop[,c(1,3:9)];
#colnames(dat_rej2) = c("Percent","RhoX","h2y","beta_prop","h2","Method","M","n2")
#dat_rej2$PropC = rep("Power", 648);
#dat_rej =rbind(dat_rej1,dat_rej2);
#dat_rej = rej_prop;
dat_rej = dat[[3]];
dat_rej$h2z=paste("",dat_rej$h2,sep="")
#dat_rej$h2y=paste("",dat_rej$h2y,sep="")
dat_rej$Power = dat_rej$rej_prop
dat_rej$Sparsity = dat_rej$beta_prop
dat_rej$sd_rej = as.numeric(as.character(dat_rej$sd_rej))
dat_rej = dat_rej[dat_rej$Method!="2-stage:AUDI",]
#dat_rej$Method[dat_rej$Method=="AUDI"] ="CoMM"
library(plyr)
dat_rej$Method=revalue(dat_rej$Method, c("AUDI"="CoMM"))
dat_rej$Method=droplevels(dat_rej$Method)
#mapvalues(dat_rej$Method, from = c("AUDI"), to = c("CoMM"))

library(ggplot2)
for (rho in c(-0.8,-0.5,-0.2,0.2,0.5,0.8)){
for (n2 in c(4000,8000)){
t1e_rej = dat_rej[dat_rej$RhoX==rho&dat_rej$n2==n2,]
#t1e_rej = dat_rej[dat_rej$RhoX==rho&dat_rej$beta_prop==beta_prop&dat_rej$n2==n2,]
#t1e_rej = rej_prop[rej_prop$RhoX==rho&rej_prop$beta_prop==beta_prop,]
#t1e_rej = rej_prop[rej_prop$RhoX==0.2&rej_prop$beta_prop==1&t1e_rej$h2==0&t1e_rej$h2y==0.4&t1e_rej$n2==5000,]
t1e_rej$h2z = factor(t1e_rej$h2z)
t1e_rej$h2y = factor(t1e_rej$h2y)
t1e_rej$Sparsity = factor(t1e_rej$Sparsity)
t1e_rej$n2 = factor(t1e_rej$n2)
t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("CoMM","2-stage:Ridge","2-stage:Enet","SKAT"))
t1e_rej$Power = as.numeric(as.character((t1e_rej$Power)))
#t1e_rej$PropC = factor(t1e_rej$PropC);
#t1e_rej$PropC = ordered(t1e_rej$PropC,levels = c("Power","NA"));

#
#rej_propu = filter(rej_prop, rank(rej_prop) > 1)
t1e_rej$h2y2 <- factor(t1e_rej$h2y, labels = c("h[C]^2==0.01", "h[C]^2==0.03", "h[C]^2==0.05", "h[C]^2==0.07", "h[C]^2==0.09"))
t1e_rej$h2z2 <- factor(t1e_rej$h2z, labels = c("h[T]^2==0", "h[T]^2==0.001", "h[T]^2==0.002", "h[T]^2==0.003"))


ggplot(t1e_rej, aes(x = Sparsity, y = Power,fill = Method))+
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=Power-sd_rej, ymax=Power+sd_rej), width=.2,
                 position=position_dodge(.9)) +
facet_grid(h2z2~h2y2,labeller = label_parsed,scales = "free_y") + 
#scale_y_continuous(limits=c(0,0.8),breaks=seq(0,1,0.2))+
#theme(panel.grid.major = element_blank(),legend.position="bottom",axis.ticks.y = element_blank())
theme(legend.position="bottom",axis.title.x = element_text(size=18),
axis.text=element_text(size=14),axis.title.y=element_blank(),
legend.title=element_text(size=20), legend.text=element_text(size=20),
strip.text.x =element_text(size=18),strip.text.y =element_text(size=18))+
geom_hline(yintercept=0.05,colour="orange",linetype="dashed")

outfile = paste("barplot2_rho_",gsub('[.]','',as.character(rho)),
"_n1_",n1,"_n2_",n2,".png",sep="")
ggsave(outfile, width = 320, height = 200, units = "mm", dpi=600)
}}


##########################################################
##########################################################
##########################################################
######type 1 error
rm(list = ls())
setwd("D:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\Simulation\\Archive\\sim")

L = 1; M = 100;
n1 = 350; #n2 = 5000;
q1 = 1; q2 = 1;

#rej_prop = array(0,dim=c(3,4,4,5,3));
#alpha_est = array(0,dim=c(3,4,4,5,2));
nrep = 1000;
rhoall = c(0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0);
n2all = c(4000,8000)

pvalall = NULL;#matrix(nrow=500*3*4*4*4*4,ncol=5)
alphaall = NULL;#matrix(nrow=500*3*4*4*4*3,ncol=5)
rej_prop = NULL;
na_prop = NULL;
for ( i in 1:3){
    for (j in 1:5){
        for (k in 1:6){
            for ( l in 1){
              for (nn in 1:2){
                rho = rhoall[i];
                h2y = h2yall[j];
                beta_prop = beta_propall[k];
                h2 = h2all[l];
                n2 = n2all[nn]

                outfile1 <- paste("sim_t1e_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                outfile2 <- paste("sim_t1e_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")

                pval <- matrix(unlist(read.table(outfile1,header=F)),ncol=5);
                alpha <- matrix(unlist(read.table(outfile2,header=F)),ncol=4);
                pval[pval==0]=1e-50
                pval[is.na(pval)] = 1;
                #rej_prop[i,j,k,l,] = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T)
                #alpha_est[i,j,k,l,] = apply(alpha, 2 , mean,na.rm=T)
                pvalall = rbind(pvalall,data.frame( as.numeric(as.character(pval)),rep(rho,nrep*5),
                        rep(h2y,nrep*5),rep(beta_prop,nrep*5), rep(h2,nrep*5),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT","2-stage:AUDI"),each=nrep),
                        rep(M,nrep*5), rep(n2,nrep*5)))
                alphaall= rbind(alphaall,data.frame( as.numeric(as.character(alpha)), rep(rho,nrep*4),
                        rep(h2y,nrep*4),rep(beta_prop,nrep*4), rep(h2,nrep*4),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","2-stage:AUDI"),each=nrep),
                        rep(M,nrep*4), rep(n2,nrep*4)))
                tmp_rej = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T);
                sd_rej =  sqrt(tmp_rej*(1-tmp_rej)/nrep)
                tmp_na_rej = apply(is.na(pval) ,2,sum) / nrep;
                rej_prop = rbind(rej_prop,data.frame(cbind(tmp_rej,tmp_na_rej,sd_rej, 
                        rep(rho,5),rep(h2y,5),rep(beta_prop,5),
                        rep(h2,5),c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT","2-stage:AUDI"), rep(M,5),rep(n2,5))))
                #na_prop = rbind(na_prop, data.frame(cbind(tmp_na_rej,rep(rho,4),rep(h2y,4),rep(beta_prop,4),
                #        rep(h2,4),c("AUDI","2-stage:Ridge","2-stage:Lasso","SKAT"), rep(M,4),rep(n2,4))))
            }
            }

        }

    }

}

pvalall = data.frame(pvalall)
colnames(pvalall) = c("Pval","RhoX","h2y","beta_prop","h2","Method","M","n2")
alphaall = data.frame(alphaall)
colnames(alphaall) = c("Alpha","RhoX","h2y","beta_prop","h2","Method","M","n2")
rej_prop = data.frame(rej_prop);
colnames(rej_prop) = c("rej_prop","na_prop","sd_rej","RhoX","h2y","beta_prop","h2","Method","M","n2")
#na_prop = data.frame(na_prop);
#colnames(na_prop) = c("na_prop","RhoX","h2y","beta_prop","h2","Method","M","n2")

dat = list();
dat[[1]] = pvalall;
dat[[2]] = alphaall;
dat[[3]] = rej_prop;

save(dat, file = "sim_t1e.RData")

##barplot
#dat_rej1 = rej_prop[,2:9]
#colnames(dat_rej1) = c("Percent","RhoX","h2y","beta_prop","h2","Method","M","n2")
#dat_rej1$PropC = rep("NA", 216);
#dat_rej2 = rej_prop[,c(1,3:9)];
#colnames(dat_rej2) = c("Percent","RhoX","h2y","beta_prop","h2","Method","M","n2")
#dat_rej2$PropC = rep("Power", 216);
#dat_rej =rbind(dat_rej1,dat_rej2);
#dat_rej$h2C=paste("h2z:",dat_rej$h2,sep=" ")
#dat_rej$h2yC=paste("h2y:",dat_rej$h2y,sep=" ")
dat_rej = rej_prop;
dat_rej$h2z=paste("",dat_rej$h2,sep="")
#dat_rej$h2y=paste("",dat_rej$h2y,sep="")
dat_rej$Power = dat_rej$rej_prop
dat_rej$Sparsity = dat_rej$beta_prop
dat_rej$sd_rej = as.numeric(as.character(dat_rej$sd_rej))

library(ggplot2)
for (rho in c(0.2,0.5,0.8)){
t1e_rej = dat_rej[dat_rej$RhoX==rho,]
t1e_rej$h2z = factor(t1e_rej$h2z)
t1e_rej$h2y = factor(t1e_rej$h2y)
t1e_rej$Sparsity = factor(t1e_rej$Sparsity)
t1e_rej$n2 = factor(t1e_rej$n2)
t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))
t1e_rej$Power = as.numeric(as.character((t1e_rej$Power)))

#t1e_rej = rej_prop[rej_prop$RhoX==rho&rej_prop$beta_prop==beta_prop,]
#t1e_rej = rej_prop[rej_prop$RhoX==0.2&rej_prop$beta_prop==1&t1e_rej$h2==0&t1e_rej$h2y==0.4&t1e_rej$n2==5000,]
#t1e_rej$h2z = factor(t1e_rej$h2z)
#t1e_rej$h2y = factor(t1e_rej$h2y)
#t1e_rej$n2 = factor(t1e_rej$n2)
#t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))
#t1e_rej$Percent = as.numeric(as.character((t1e_rej$Percent)))
#t1e_rej$PropC = factor(t1e_rej$PropC);
#t1e_rej$PropC = ordered(t1e_rej$PropC,levels = c("Power","NA"));

#rej_propu = filter(rej_prop, rank(rej_prop) > 1)
ggplot(t1e_rej, aes(x = Sparsity, y = Power,fill = Method))+
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=Power-sd_rej, ymax=Power+sd_rej), width=.2,
                 position=position_dodge(.9)) +
facet_grid(n2~h2y,labeller = label_both,scales = "free_y") + 
scale_y_continuous(breaks=seq(0,0.2,0.05),limits=c(0,0.2))+
#theme(panel.grid.major = element_blank(),legend.position="bottom",axis.ticks.y = element_blank())
theme(legend.position="bottom",axis.title.x =element_text(size=18),
axis.text=element_text(size=18),axis.title.y=element_blank(),
legend.title=element_text(size=20), legend.text=element_text(size=20),
strip.text.x =element_text(size=18),strip.text.y =element_text(size=18))+
geom_hline(yintercept=0.05,colour="orange",linetype="dashed")

outfile = paste("barplot_t1e_rho_",gsub('[.]','',as.character(rho)),".png",sep="")
ggsave(outfile, width = 320, height = 240, units = "mm", dpi=600)
}

##########################################################
##########################################################
##########################################################
##########test
rm(list = ls())
setwd("G:\\Data\\GoogleDrive\\Work\\Multiple_Trait\\TWAS\\AUDI\\Simulation")


L = 1; M = 100;
n1 = 350; #n2 = 5000;
q1 = 1; q2 = 1;
nrep = 50;

#rej_prop = array(0,dim=c(3,4,4,5,3));
#alpha_est = array(0,dim=c(3,4,4,5,2));
rhoall = c(0.2,0.5,0.8);
h2yall = c(0.1,0.3, 0.5)
beta_propall = c(0.01, 0.05, 0.1, 0.5)
h2all = c(0.0, 0.001, 0.005, 0.01);
n2all = c(10000,20000)

pvalall = NULL;#matrix(nrow=500*3*4*4*4*4,ncol=5)
alphaall = NULL;#matrix(nrow=500*3*4*4*4*3,ncol=5)
rej_prop = NULL;
for ( i in 1){
    for (j in 2){
        for (k in 1:3){
            for ( l in 2){
              for (nn in 1){
                rho = rhoall[i];
                h2y = h2yall[j];
                beta_prop = beta_propall[k];
                h2 = h2all[l];
                n2 = n2all[nn]

                outfile1 <- paste("sim_AUDI_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                outfile2 <- paste("sim_AUDI_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")

                pval <- matrix(unlist(read.table(outfile1,header=F)),ncol=4);
                alpha <- matrix(unlist(read.table(outfile2,header=F)),ncol=3);
                #rej_prop[i,j,k,l,] = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T)
                #alpha_est[i,j,k,l,] = apply(alpha, 2 , mean,na.rm=T)
                pvalall = rbind(pvalall,data.frame( as.numeric(as.character(pval)),rep(rho,nrep*4),
                        rep(h2y,nrep*4),rep(beta_prop,nrep*4), rep(h2,nrep*4),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"),each=nrep),
                        rep(M,nrep*4), rep(n2,nrep*4)))
                alphaall= rbind(alphaall,data.frame( as.numeric(as.character(alpha)), rep(rho,nrep*3),
                        rep(h2y,nrep*3),rep(beta_prop,nrep*3), rep(h2,nrep*3),
                        rep(c("AUDI","2-stage:Ridge","2-stage:Enet"),each=nrep),
                        rep(M,nrep*3), rep(n2,nrep*3)))
                tmp_rej = apply(pval < 0.05 ,2,sum,na.rm=T) / apply(!is.na(pval),2,sum,na.rm=T);
                rej_prop = rbind(rej_prop,data.frame(cbind(tmp_rej,rep(rho,4),rep(h2y,4),rep(beta_prop,4),
                        rep(h2,4),c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"), rep(M,4),rep(n2,4))))
            }
            }

        }

    }

}

pvalall = data.frame(pvalall)
colnames(pvalall) = c("Pval","RhoX","h2y","beta_prop","h2","Method","M","n2")
alphaall = data.frame(alphaall)
colnames(alphaall) = c("Alpha","RhoX","h2y","beta_prop","h2","Method","M","n2")

rej_prop = data.frame(rej_prop);
colnames(rej_prop) = c("rej_prop","RhoX","h2y","beta_prop","h2","Method","M","n2")

alphaall$beta_prop = factor(alphaall$beta_prop)
alphaall$Method <- ordered(alphaall$Method, levels = c("AUDI","2-stage:Ridge","2-stage:Enet","SKAT"))

library(ggplot2)
ggplot(pvalall, aes(x =h2  , y = -log10(Pval),fill =Method ))+ geom_boxplot()+
geom_hline(yintercept=0.3,colour="orange",linetype="dashed") +
facet_grid(beta_prop~h2y,labeller = label_both,scales = "free_y") +
theme(legend.position="bottom")


