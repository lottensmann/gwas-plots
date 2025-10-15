
bgen_file_prefix="imputed_chr" #prefix of bgen file (before chromosome number)
bgen_file_suffix="_SNPids.bgen" #suffix of bgen file (after chromosome number)
bucket="gs://genotypes/" #gcloud bucket for bgen files
gt_path="lottensm/gt/" #location of bgen files
vep_path="lottensm/vep/" #location of output files
output_name="snps" #name of output files

#read in file with lead variants
#required columns chr, variant (chr:bp:ref:alt), variant_id ("chr"chr_bp_ref_alt), rsid
leadvars=read.delim("leadvars.txt",header=T,sep="\t",stringsAsFactors = F)

for(i in unique(leadvars$chr)){
  leadvars_chr=leadvars[which(leadvars$chr==i),]
  #download bgen file from gcloud bucket
  if(!file.exists(paste0(gt_path,bgen_file_prefix,i,bgen_file_suffix))){
    system(paste0("gsutil cp ",bucket,bgen_file_prefix,i,bgen_file_suffix," ",gt_path))
  }
  #create file with variant ids as input for qctool
  snps=unique(leadvars_chr$variant_id) 
  write.table(snps,file=paste0(vep_path,output_name,"_chr",i,".txt"), row.names=F, col.names=F, sep='\t', quote=F)
  #create vcf file with qctool
  system(paste0("qctool -g ",gt_path,bgen_file_prefix,i,bgen_file_suffix," -incl-snpids ",vep_path,output_name,"_chr",i,".txt -og ",vep_path,output_name,"_chr",i,".vcf"))
  #subset for the first 5 columns
  system(paste("cat ",vep_path,output_name,"_chr",i,".vcf | vcf-subset -c CHROM,POS,ID,REF,ALT -t ref,SNPs,indels,MNPs,other -u > ",vep_path,output_name,"_chr_",i,"_subcols.vcf",sep=""))
  #delete bgen file 
  system(paste0("rm ",gt_path,bgen_file_prefix,i,bgen_file_suffix))
}

#concat vcf files from all chromosomes
system(paste("vcf-concat ",vep_path,output_name,"_chr_*_subcols.vcf > ",vep_path,output_name,"_all_chr_subcols.vcf",sep=""))
#check number of included variants
system("wc -l ",vep_path,output_name,"_all_chr_subcols.vcf") 

#get nearest genes with vep, default distance 5000
distance=300000

system(paste0("vep --force_overwrite --database -i ",vep_path,output_name,"_all_chr_subcols.vcf -o ",
              vep_path,output_name,"_vep_genes_protein.txt --distance ",distance,",",distance," --uniprot --symbol"))
vep=read.delim(paste0(vep_path,output_name,"_vep_genes_protein.txt"),header=F,stringsAsFactors=F, sep = '\t',comment.char = "#") 

vep$distance=0
vep$symbol=""
for(i in 1:dim(vep)[1]){
  #get distance from output file
  split_dist=unlist(strsplit(vep[i,]$V14,"DISTANCE=")) 
  if(length(split_dist)>1){
    #use first listed distance
    split_dist2=unlist(strsplit(split_dist[2],";")) 
    #add to data frame
    vep[i,]$distance=as.numeric(split_dist2[1]) 
  }
  #get uniprot symbol from output file
  split_symbol=unlist(strsplit(vep[i,]$V14,"SYMBOL=")) 
  if(length(split_symbol)>1){
    #use first listed symbol
    split_symbol2=unlist(strsplit(split_symbol[2],";")) 
    vep[i,]$symbol=split_symbol2[1] 
  }
}
#save file with distances
write.table(vep,file=paste0(vep_path,output_name,"_vep_distance.txt"), row.names=F, col.names=T, sep='\t', quote=F)


#add all nearby genes sorted by distance to data frame of variants
leadvars$nearest_genes="" #uniprot gene symbol
leadvars$nearest_ensg="" #ensg symbol
leadvars$nearest_distance="" #distance

for(i in 1:dim(leadvars)[1]){
  var=leadvars[i,]$variant_id
  vep_snp=vep[which(vep$V1==var),] #all genes for the variant
  vep_snp=vep_snp[order(vep_snp$distance),] #order by distance
  vep_snp=vep_snp[!duplicated(vep_snp$V4),]
  #collape distances and genes and add to data frame
  leadvars[i,]$nearest_genes=paste(vep_snp$symbol,collapse = ";")
  leadvars[i,]$nearest_ensg=paste(vep_snp$V4,collapse = ";")
  leadvars[i,]$nearest_distance=paste(vep_snp$distance,collapse = ";")
}
write.table(leadvars,file=paste0(vep_path,output_name,"_vep_nearest_genes.txt"), row.names=F, col.names=T, sep='\t', quote=F)
