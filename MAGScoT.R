### MAGScoT - Bin Scoring and Refinment Tool for Metagenome Assblies

library("optparse")

option_list = list(
		make_option(c("-i", "--input"), type="character", default=NULL,
							help="Tab-separated input file with three columns: bin, contig, set; no header!", metavar="character"),
		make_option(c("--hmm"), type="character", default=NULL,
							help="Tab-separated input file with marker mapping; three columns: gene id, marker, e-value. \n Gene IDs must represent contig IDs after removal of _[0-9] at the end of the name (prodigal default)", metavar="character"),
  	make_option(c("-p", "--profile"), type="character", default="default",
              help="Profile used for scoring, all derived from GTDB release 207. default[=bac120+ar53], ar53, bac120. [default]", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="MAGScoT",
              help="output file name base [default=%default]", metavar="character"),
		make_option(c("-a", "--score_a"), type="double", default=1, help="Scoring parameter a [default=%default]"),
		make_option(c("-b", "--score_b"), type="double", default=0.5, help="Scoring parameter b [default=%default]"),
		make_option(c("-c", "--score_c"), type="double", default=0.5, help="Scoring parameter c [default=%default]"),
		make_option(c("-t", "--threshold"), type="double", default=0.5, help="Scoring minimum completeness threshold [default=%default]"),
		make_option(c("--score_only"), action="store_true", dest="score_only", help="Only do scoring, no refinement [false]"),
		make_option(c("--skip_merge_bins"), action="store_true", dest="skip_merge_bins", help="Skip bin merging [false]"),
		make_option(c("-m", "--min_markers"), type="double", default=25, help="Minimum number of unique markers in bins to be considered as seed for bin merging [default=%default]"),
		make_option(c("-s", "--min_sharing"), type="double", default=0.8, help="Minimum percentage of shared markers for bin sharing. [default=%default]"),
		make_option(c("-n", "--n_iterations"), type="double", default=2, help="Number of merging iterations to perform. [default=%default]")

);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat("Loading packages...\n")
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(funr))
### READ IN BINNiNG INFO
### FORMAT: 3 columns => bin, contig, set

contig_to_bin<-read_delim(file=opt$input, delim="\t", col_names=c("bin","contig","set"), show_col_types = FALSE)
contig_to_bin=contig_to_bin %>% filter(!is.na(set))
if(ncol(contig_to_bin) != 3){
	cat("Contig to bin input file should have three tab-separated colums without a header: Bin ID, Contig ID, Binning Tool\n")
	cat("Your file ",opt$input, " has ",ncol(contig_to_bin)," tab-separated columns.\n")
	cat("Please check your input files and run MAGScoT again.\nAborting execution!\n")
	quit(save="no")
}

if(length(unique(contig_to_bin$bin)) > length(unique(contig_to_bin$contig))){
	cat("Contig to bin input file should have three tab-separated colums without a header: Bin ID, Contig ID, Binning Tool\n")
	cat("The number of distinct items in the column 'Contig ID' is smaller than in the column 'Bin ID'.\n Is your input correctly formatted?\n")
	cat("Please check your input files and run MAGScoT again.\nAborting execution!\n")
	quit(save="no")
}

num_binsets = length(unique(contig_to_bin$set))

if(num_binsets==1){
	cat("Only output from one binning alogorithm provided. Skipping merging and refinement. Scoring bins only.\n")
	opt$score_only = T
	opt$n_iter = 0
}

if(opt$n_iter >= num_binsets & is.null(opt$skip_merge_bins) == T){
	cat("Number of bin-merging iterations cannot be equal or larger than number or binning sofwares used.\n")
	cat("Setting iterations to:", num_binsets-1,"\n")
	opt$n_iter = num_binsets -1
}

id=opt$out

### Parameter setup a=1,b=0.5,c=0.5 is DAS Tool default; metawrap is much more stringent on contamination (b=5)
a=opt$score_a
b=opt$score_b
c=opt$score_c
cutoff=opt$threshold

### Parameters for bin merging
min_markers = opt$min_markers
min_sharing = opt$min_sharing

#### Read markers file

if(opt$profile=="default"){
	markers=read.table(file.path(funr::get_script_path(),"profiles","gtdb_rel207_default_markers.tsv"), stringsAsFactors=F, head=T)
}else if(opt$profile=="bac120"){
	markers=read.table(file.path(funr::get_script_path(),"profiles","gtdb_rel207_bac120_markers.tsv"), stringsAsFactors=F, head=T)
}else if(opt$profile=="ar53"){
	markers=read.table(file.path(funr::get_script_path(),"profiles","gtdb_rel207_ar53_markers.tsv"), stringsAsFactors=F, head=T)
}else{
	print_help(opt_parser)
  stop("No valid profile selected.n", call.=FALSE)
}
markers = markers %>% mutate(marker=gsub("[.]HMM|[.]hmm","", marker))


### READ IN Singe-copy gene annotation
### only best annotation per gene is kept => sort by X3 (evalue) and remove duplicates
### if the same gene occurs twice in a contig it is only counted once
hmm<-read_delim(opt$hmm, delim="\t",col_names=FALSE, show_col_types = FALSE)
if(ncol(hmm) != 3){
	cat("HMM input file should have three tab-separated colums without a header: Protein ID, Marker ID, e-value\n")
	cat("Your file ",opt$hmm, " has ",ncol(hmm)," tab-separated columns.\n")
	cat("Please check your input files and run MAGScoT again.\nAborting execution!\n")
	quit(save="no")
}

gene_to_contig=hmm %>% arrange(X3) %>%
	filter(!duplicated(X1)) %>%
	dplyr::select(X1, X2) %>%
	mutate(contig=gsub("_[0-9]+$","",X1), gene=X2) %>%
	dplyr::select(contig, gene) %>%
	distinct()

####
if(is.null(opt$skip_merge_bins) == F | is.null(opt$score_only) == F | opt$n_iter == 0){
	cat("Skipping bin merging...\n")
}else{
	for(it in 1:opt$n_iter){
		cat("Bisco iteration: ", it, "\n")
		if(it==1){
			contig_to_bin_bisco_it = contig_to_bin
			contig_to_bin_bisco = contig_to_bin
		}
		gene_to_bin = left_join(contig_to_bin_bisco_it %>% filter(contig %in% gene_to_contig$contig), gene_to_contig, by="contig")
		merge_cand = gene_to_bin %>% distinct_at(.vars=c("bin","gene"), .keep_all=T) %>% group_by(bin) %>% summarise(n_markers = n()) %>%
			left_join( gene_to_bin %>% group_by(bin) %>% summarise(n_markers_tot = n()), by="bin") %>% filter(n_markers >= min_markers, n_markers_tot <= 150) %>% pull(bin) %>% unique()
		if(length(merge_cand)==0) break
		binmatch=sapply(seq_along(merge_cand), function(thisbin_id){
			cat("Processing candidate bin ",thisbin_id, " of ", length(merge_cand),"\r")
			thisbin = merge_cand[thisbin_id]
			thisbin_contigs = contig_to_bin_bisco_it %>% filter(bin == thisbin, contig %in% gene_to_contig$contig) %>% pull(contig)
			thisbin_genes = gene_to_contig %>% filter(contig %in% thisbin_contigs)
			if(it>1){thisbin_blacklist = bisco_out %>% filter(bin==thisbin) %>% select(bin_a, bin_b) %>% unlist %>% as.vector()}else{thisbin_blacklist=c()}
			return(contig_to_bin %>% filter(contig %in% thisbin_genes$contig) %>% left_join(thisbin_genes, by="contig") %>%
				distinct_at(.vars=c("bin","gene"), .keep_all=T) %>% group_by(bin) %>% summarize(shared_markers=n()) %>%
				mutate(shared_markers_rel = shared_markers/length(unique(thisbin_genes$gene))) %>% filter(bin != thisbin, !bin %in% thisbin_blacklist, shared_markers_rel >= min_sharing) %>%
				mutate(seed=thisbin, bin_a=ifelse(seed>bin, seed, bin), bin_b=ifelse(seed>bin, bin, seed)) %>% select(seed, bin_a, bin_b, shared_markers, shared_markers_rel))
		}, simplify=F)
		bisco_out = binmatch %>% do.call("rbind",.) %>% distinct_at(.vars=c("bin_a","bin_b"), .keep_all=T) %>% filter(shared_markers >= min_markers) %>% mutate(bin=paste0(paste0("MAGScoT_",it,"_"), seq_along(shared_markers)))
		contig_to_bin_bisco_it = bisco_out %>%	apply(., 1, function(x) contig_to_bin %>% filter(bin %in% c(x[2], x[3])) %>%
			mutate(bin=x[6], set=paste0("MAGScoT_",it)) %>% distinct_at(.vars=c("bin","contig"), .keep_all=T)) %>% do.call("rbind", .) %>% data.frame(row.names=NULL)
		contig_to_bin_bisco = rbind(contig_to_bin_bisco, contig_to_bin_bisco_it)
		if(it==1){bisco_out_all = bisco_out}else{bisco_out_all = rbind(bisco_out_all, bisco_out)}
		if(nrow(contig_to_bin_bisco_it)==0){cat("Exiting bin merging, no (more) overlaps found\n"); break}
		cat("\n")
	}
	contig_to_bin = contig_to_bin_bisco
}

####
contig_to_bin.remain<-contig_to_bin
contig_to_bin.out<-contig_to_bin[0,]

### Initial values
i=1
allmax=100
recalc = contig_to_bin.remain$bin %>% unique

while(nrow(contig_to_bin.remain)>1){
cat("Refining bins, iteration:", i, "\n")
cat("Extracting SCG information for bins...\n")

### at first run (exist("zz") == F) all gene_to_bin statistics (zz) are calculated
### at later iterations, only bins with changed contents are re-calculated; this saves a lot of time

if(length(recalc) > 0){
	gene_to_bin = left_join(contig_to_bin.remain %>% filter(contig %in% gene_to_contig$contig), gene_to_contig, by="contig")
	if(exists("zz")){
		zz = zz %>% filter(bin %in% scoreframe$bin, !bin %in% recalc)
		zz = rbind(zz, gene_to_bin %>% filter(bin %in% recalc) %>% select(bin, gene) %>% group_by(bin,gene) %>% summarize(count=n(), .groups="drop"))
	}else{
		zz=gene_to_bin %>% select(bin, gene) %>% group_by(bin,gene) %>% summarize(count=n(), .groups="drop")
	}

} else {
zz = zz %>% filter(bin %in% scoreframe$bin)
}

### stat calculation for marker sets; also here a more flexible approch would be great, e.g. iterating through a set of markers
df_list = sapply(unique(markers$set), function(this_set){
	SET_MARKERS = markers %>% filter(set == this_set) %>% pull(marker)
	zz %>% filter(gene %in% SET_MARKERS) %>%
	group_by(bin) %>% summarize(uniqueSCGs = n(), multipleSCGs = sum(count>1), sumSCGs = sum(count)) %>%
	mutate(Completeness = uniqueSCGs / length(SET_MARKERS), additionalSCGs = sumSCGs - uniqueSCGs - multipleSCGs, Contamination = b*(multipleSCGs / uniqueSCGs) + c*(additionalSCGs/length(SET_MARKERS)), score = a*Completeness - Contamination)}, simplify=F)

if(i==1){for(set in names(df_list)){write.table(df_list[set], paste0(id,".",set,".out"), sep="\t", row.names=F)}}

if(i==1){
	scoreframe = data.frame(bin=unique(contig_to_bin$bin), stringsAsFactors=F)
}else{
	scoreframe = data.frame(bin=unique(zz$bin), stringsAsFactors=F)
}
for(set in names(df_list)){
	scoreframe = scoreframe %>% left_join(df_list[[set]] %>% rename_with(~paste0(colnames(df_list[[set]])[-1],".", set), all_of(colnames(df_list[[set]])[-1])), by="bin")
	scoreframe[is.na(scoreframe)] = 0
}

scoreframe$max = scoreframe %>% select(contains("score.")) %>% apply(., 1, max, na.rm=T)
scoreframe$max_complete = scoreframe %>% select(contains("Completeness.")) %>% apply(., 1, max, na.rm=T)

scoreframe = scoreframe %>% arrange(-max, set) %>% left_join(contig_to_bin %>% dplyr::select(bin,set) %>% distinct, by="bin")

if(i==1){
	write.table(scoreframe, paste0(id,".scores.out"), sep="\t", row.names=F)
	cat("Scores for all initial bins written to:" , paste0(id,".scores.out"), "\n")
	if(!is.null(opt$score_only)){
		quit(save="no")
	}
}

if(nrow(scoreframe)==0 | scoreframe[1,"max"] < cutoff){
cat("No bin surpasses the cutoff of:" , cutoff, "\n")
break
}

### bins where neither bac not arc markers are above the cutoff (50% by DAS tool default) are removed as they will never reach the scoring cutoff
scoreframe = scoreframe %>% filter(max_complete >= cutoff)

### as long as consecutive best hits are from the same set, no high-scoring bins are affected by definition (each contig only occurs once per set)
### this iteration can mean big speedup

winnerset=scoreframe$set[1]

while(nrow(scoreframe)>0 & scoreframe$set[1] == winnerset){
winnerbin = scoreframe$bin[1]

if(exists("scoreframe.out")==F){scoreframe.out<-head(scoreframe, 1)}else{scoreframe.out<-rbind(scoreframe.out,head(scoreframe,1))}

cat("The highest scoring bin is:", winnerbin,"\n")
### add contigs for winning bin to output table
contig_to_bin.out <- rbind(contig_to_bin.out,contig_to_bin.remain %>% filter(bin==winnerbin))

#print(tail(scoreframe.out,1))
cat("Remaining:", length(unique(contig_to_bin.remain$contig)), "contigs and", nrow(scoreframe),"candidate bins\n")
scoreframe=scoreframe[-1,]
if(substr(winnerset, 1, 7)=="MAGScoT") break
}

### assess which bins were affected by the winning bin(s); store these for re-scoring in the next iteration
recalc = contig_to_bin.remain %>% filter(contig %in% contig_to_bin.out$contig) %>% pull(bin)

### remove bins below treshold (already removed in df); remove contigs which have been assigned to a refined bin
contig_to_bin.remain <- contig_to_bin.remain %>% filter(bin %in% scoreframe$bin, !contig %in% contig_to_bin.out$contig)

i=i+1
}

if(exists("scoreframe.out")){
### rename refined bins
	rownames(scoreframe.out)<-paste0(id,"_cleanbin_",formatC(seq_along(scoreframe.out$bin), width = 6, format = "d", flag = "0"))
	contig_to_bin.out$binnew<-rownames(scoreframe.out)[match(contig_to_bin.out$bin, scoreframe.out$bin)]

###
	cat("Refinement lead to a total of", nrow(scoreframe.out)," bins with a score >=", cutoff,"\n")
	stats_outfile = paste0(id,".refined.out")
	binning_outfile = paste0(id,".refined.contig_to_bin.out")

	cat("Refinement stats are written to:", stats_outfile,"\n")
	cat("Contig-to-refined-bin mapping is written to:", binning_outfile, "\n")

	write.table(scoreframe.out,stats_outfile,sep="\t",quote=F)
	write.table(contig_to_bin.out %>% dplyr::select(binnew,contig),binning_outfile,sep="\t",quote=F, row.names=F)
} else {
	cat("No bins with score >=", cutoff, "were found in the dataset.\n")
}
