library(ape)
library(phangorn)
set.seed(0)

#Load the phylogenetic tree with clusters
t=read.tree('clusters-covid.nwk')
t$root.time=2020.462
s=unlist(strsplit(t$tip.label,'-'))
membership=as.numeric(s[seq(2,length(s),2)])

#Plot the tree with clusters
dates=t$root.time+unname(dist.nodes(t)[Ntip(t)+1,])
parents=c();parents[t$edge[,2]]=t$edge[,1]
wi=c();wi[t$edge[,2]]=1:Nedge(t)
col=rep(1,Nedge(t))
for (i in 1:length(membership)) {k=i;while (dates[k]>2021) {
  col[wi[k]]=membership[i]+1
  k=parents[k]}
}
pdf('clusters-covid.pdf',10,10)
palette(c('#000000',rep(rainbow(10),100)))
plot(t,show.tip.label = F,tip.color=membership,cex=0.2,align.tip.label = F,edge.color=col)
axisPhylo(1,backward = F)
dev.off()

#Generate the WGS data
rootseq=read.dna('wuhan.fasta',format='fasta',as.character = T)
d=simSeq(t,rate=1e-3,rootseq=rootseq,l=length(rootseq))

system('mkdir Train;mkdir Test')
for (i in 1:max(membership)) system(sprintf('mkdir Train/%d;mkdir Test/%d',i,i))

for (i in 1:Ntip(t)) {
  w=which(membership==membership[i])
  if (i<=w[floor(length(w)/2)]) type='Train' else type='Test'
  write.phyDat(d[i],sprintf('%s/%d/%s.fasta',type,membership[i],t$tip.label[i]),format='fasta')
}
