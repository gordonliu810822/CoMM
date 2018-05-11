rm(list = ls())

L = 1; M = 100;
n1 = 250; #n2 = 5000;
q1 = 1; q2 = 1;

nrep = 500;
rhoall = c(-0.8,-0.5,-0.2,0.2,0.5,0.8);
h2yall = c(0.01,0.03, 0.05, 0.07,0.09)
beta_propall = c(0.1,0.2,0.3,0.4,0.5,1)
h2all = c(0, 0.001, 0.002, 0.003);
n2all = c(4000,8000)

pvalall = NULL;
alphaall = NULL;
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
                    outfile1 <- paste("sim_t1e_CoMM_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                    outfile2 <- paste("sim_t1e_CoMM_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
		}
                if (h2 != 0) {
                    outfile1 <- paste("sim_CoMM_pval_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
                         q2,"_",h2y,"_",beta_prop,"_",h2,".txt",sep="")
                    outfile2 <- paste("sim_CoMM_alpha_",n1,"_",n2,"_",L,"_",M,"_",rho,"_",q1,"_",
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

dat = list();
dat[[1]] = pvalall;
dat[[2]] = alphaall;
dat[[3]] = rej_prop;

save(dat, file = "sim3_n1_250.RData")

dat_rej = dat[[3]];
dat_rej$h2z=paste("",dat_rej$h2,sep="")
dat_rej$Power = dat_rej$rej_prop
dat_rej$Sparsity = dat_rej$beta_prop
dat_rej$sd_rej = as.numeric(as.character(dat_rej$sd_rej))
dat_rej = dat_rej[dat_rej$Method!="2-stage:AUDI",]
library(plyr)
dat_rej$Method=revalue(dat_rej$Method, c("AUDI"="CoMM"))
dat_rej$Method=revalue(dat_rej$Method, c("2-stage:Ridge"="PrediXcan:Ridge"))
dat_rej$Method=revalue(dat_rej$Method, c("2-stage:Enet"="PrediXcan:Enet"))
dat_rej$Method=droplevels(dat_rej$Method)


library(ggplot2)
for (rho in c(-0.8,-0.5,-0.2,0.2,0.5,0.8)){
for (n2 in c(4000,8000)){
t1e_rej = dat_rej[dat_rej$RhoX==rho&dat_rej$n2==n2,]
t1e_rej$h2z = factor(t1e_rej$h2z)
t1e_rej$h2y = factor(t1e_rej$h2y)
t1e_rej$Sparsity = factor(t1e_rej$Sparsity)
t1e_rej$n2 = factor(t1e_rej$n2)
t1e_rej$Method <- ordered(t1e_rej$Method, levels = c("CoMM","PrediXcan:Ridge","PrediXcan:Enet","SKAT"))
t1e_rej$Power = as.numeric(as.character((t1e_rej$Power)))

t1e_rej$h2y2 <- factor(t1e_rej$h2y, labels = c("h[C]^2==0.01", "h[C]^2==0.03", "h[C]^2==0.05", "h[C]^2==0.07", "h[C]^2==0.09"))
t1e_rej$h2z2 <- factor(t1e_rej$h2z, labels = c("h[T]^2==0", "h[T]^2==0.001", "h[T]^2==0.002", "h[T]^2==0.003"))


ggplot(t1e_rej, aes(x = Sparsity, y = Power,fill = Method))+
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=Power-sd_rej, ymax=Power+sd_rej), width=.2,
                 position=position_dodge(.9)) +
facet_grid(h2z2~h2y2,labeller = label_parsed,scales = "free_y") + 
theme(legend.position="bottom",axis.title.x = element_text(size=18),
axis.text=element_text(size=14),axis.title.y=element_blank(),
legend.title=element_text(size=20), legend.text=element_text(size=20),
strip.text.x =element_text(size=18),strip.text.y =element_text(size=18))+
geom_hline(yintercept=0.05,colour="orange",linetype="dashed")

outfile = paste("barplot2_rho_",gsub('[.]','',as.character(rho)),
"_n1_",n1,"_n2_",n2,".png",sep="")
ggsave(outfile, width = 320, height = 200, units = "mm", dpi=600)
}}



