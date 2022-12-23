# input, myabun: abundance table, n samples x m genomens
# input, mydata: metadata table, n samples x s variables

### Beta Diveristy Analysis ###
mybray=vegdist(myabun,method = "bray")
mypcoa=pcoa(mybray)
a=diag(t(mypcoa$vectors)%*%mypcoa$vectors)
explain=a/sum(a)*100
explain=round(explain,2)

### PERMANOVA test ###
set.seed(315) 
adonis2(as.matrix(mybray)~Type,mydata) #single factor
set.seed(315)
adonis2(as.matrix(mybray)~Age,mydata) # single factor
set.seed(315)
adonis2(as.matrix(mybray)~Type+Age,mydata,by="margin") # margin

### RDA analysis ###
myabun_hellinger=decostand(myabun,method = "hellinger")
decorana(myabun_hellinger)
tmpmod=rda(myabun_hellinger~Type2,mydata)
anova(tmpmod)
RsquareAdj(tmpmod)
tmpgoodness=goodness(tmpmod)

### GMI calculation ###
name1=rownames(RDAselect_33SigGenomes)[which(RDAselect_33SigGenomes$Group=="Guild1")]
name2=rownames(RDAselect_33SigGenomes)[which(RDAselect_33SigGenomes$Group=="Guild2")]
a1=apply(myabun[,name1],1,function(x){sum(x)/length(which(x>0))})
a2=apply(myabun[,name2],1,function(x){sum(x)/length(which(x>0))})
mydata$index=a1-a2

### AUROC ###
g=roc(predictor=mydata$GMI,response =mydata$Type,auc = TRUE)
ci95=plot.roc(x=mydata$Type,predictor = mydata$GMI,ci=TRUE,print.auc=TRUE)

### AUPRC ###
sscurves=pr.curve(scores.class0 =mydata$Index[mydata$Type=="group1"],
                  scores.class1 =mydata$Index[mydata$Type2=="group2"],curve=T)
AUPRC_draw=data.frame(sscurves$curve)
colnames(AUPRC_draw)=c("Recall","Precision","cutoff")