###### Arguments ######
# 1: trait file with columns FID, IID, trait
# 2: sibling relatedness table with columns: sib1_ID, sib2_ID, relatedness
# 3: output directory. outputs files trait_herit.txt and trait_n.txt to this directory
args=commandArgs(trailingOnly=T)
trait=args[1]
rel_table=args[2]
outdir = args[3]
sib_rel=read.table(rel_table,header=T, colClasses=c('character','character','numeric'))
phen=read.table(trait)
y=phen[,3]
names(y)=phen[,1]

get_phen=function(sib,y){if (sib%in%names(y)){return(y[[sib]])} else {return(NA)}}

y_sib=apply(sib_rel[,1:2],c(1,2),get_phen,y)
unique_names=unique(c(sib_rel[,1],sib_rel[,2]))
y_unique=y[names(y)%in%unique_names]
y_sib=(y_sib-mean(y_unique,na.rm=T))/sd(y_unique,na.rm=T)
y_sib=data.frame(sib_rel,y_sib,(y_sib[,1]-y_sib[,2])^2)
dimnames(y_sib)[[2]]=c('sib1','sib2','relatedness','y1','y2','y1y2')
write.table(sum(!is.na(y_sib$y1y2)),paste(outdir,trait,'_n.txt',sep=''),quote=F,row.names=F,col.names=F)

r=lm(y1y2~relatedness,data=y_sib,na.action='na.exclude')
rcoef=summary(r)$coefficients
rcoef[2,1]=-0.5*rcoef[2,1]
rcoef[2,2]=0.5*rcoef[2,2]
write.table(summary(r)$coefficients,paste(outdir,trait,'_herit.txt',sep=''),quote=F)