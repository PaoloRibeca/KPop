library(ape)
library(phangorn)

#Load the phylogenetic tree with clusters
t=read.tree('clusters-tb.nwk')
s=unlist(strsplit(t$tip.label,'-'))
membership=as.numeric(s[seq(2,length(s),2)])

#Plot the tree with clusters
pdf('clusters-tb.pdf',10,10)
grDevices::palette(grDevices::rainbow(max(membership)))
t2=t
t2$tip.label=rep('-',Ntip(t2))
plot(t2,show.tip.label = T,tip.color=membership,align.tip.label=T)
axisPhylo(1,backward = F)
dev.off()

#Generate simulated WGS data
set.seed(0)
rootseq=read.dna('refTB.fasta',format='fasta',as.character = T)
pos=which(runif(length(rootseq))<0.1)
rootseq2=rootseq[pos]
d=simSeq(t,rate=1.1e-6,rootseq=rootseq2,l=length(rootseq2))
snps=as.character(d)
snps=toupper(snps)
ref=toupper(rootseq)
for (i in 1:nrow(snps)) {
  f=file(sprintf('tmp%d.fasta',i),'w')
  print(i)
  ref[pos]=snps[i,]
  writeLines(sprintf('> %s',rownames(snps)[i]),f)
  writeLines(paste(ref,collapse = ''),f)
  close(f)
}

#Generate simulated NGS data
n=nrow(snps)
system('rm -rf Test Train Test2 Train2;mkdir Train;mkdir Test;mkdir Train2;mkdir Test2')
for (i in 1:10) system(sprintf('mkdir Train/%d;mkdir Test/%d;mkdir Train2/%d;mkdir Test2/%d',i,i,i,i))
parallel::mclapply(1:n,function (i){
  print(i)
  l=strsplit(rownames(snps)[i],'-')
  b=as.integer(l[[1]][2])
  system(sprintf('./art_illumina -rs %d -ss HS25 -i tmp%d.fasta -p -l 150 -f 20 -m 200 -s 10 -o tmp%d',i,i,i))#generate reads
  system(sprintf('rm tmp%d1.aln tmp%d2.aln',i,i))
  if (runif(1)<0.5) di='Train' else di='Test'
  system(sprintf('mv tmp%d1.fq %s/%d/%s_1.fastq;mv tmp%d2.fq %s/%d/%s_2.fastq',i,di,b,rownames(snps)[i],i,di,b,rownames(snps)[i]))
  system(sprintf('mv tmp%d.fasta %s2/%d/%s.fasta',i,di,b,rownames(snps)[i]))
}, mc.cores=16)

