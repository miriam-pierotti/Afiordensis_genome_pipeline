library(bambu)

annotation = prepareAnnotations("/nfs/scratch/pierotmi/projects/black_coral/ref_genome/4/output/braker/braker.gtf") #output from Braker
head(annotation)

se <- bambu(
    reads = "/nfs/scratch/pierotmi/projects/black_coral/rna/output/reads_to_ref/reads_1kbp_q10_to_refgenome_sorted.bam", #Reads aligned to reference genome
    annotations = annotation,
    genome = "/nfs/scratch/pierotmi/projects/black_coral/ref_genome/4/output/final_assembly/AntFio1.fasta", #Reference genome
    discovery = TRUE,
    quant = TRUE,
    verbose = TRUE,
    opt.discovery = list(fitReadClassModel = TRUE), #if TRUE, bambu perform the trainig, if FALSE, use the default model
    NDR = 0.4 # set the Novel Discovery Rate (NDR). 1 means that all novel transcripts will be retained in the output.
)

saveRDS(se, "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu_training/bambu_se.rds")
rd <- as.data.frame(rowData(se))
rd_simple <- rd[, sapply(rd, function(x) !is.list(x))]
write.csv(rd_simple, "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu_training/transcript_metadata.csv")

# NDR summary statistics (novel transcripts only)
ndr <- rd_simple$NDR[!is.na(rd_simple$NDR)]
ndr_summary <- c(
    "=== NDR Summary ===",
    paste("Total novel transcripts:", length(ndr)),
    paste("Mean NDR:  ", round(mean(ndr), 3)),
    paste("Median NDR:", round(median(ndr), 3)),
    paste("NDR < 0.1: ", sum(ndr < 0.1)),
    paste("NDR < 0.2: ", sum(ndr < 0.2)),
    paste("NDR < 0.3: ", sum(ndr < 0.3)),
    paste("NDR < 0.4: ", sum(ndr < 0.4)),
    paste("NDR < 0.5: ", sum(ndr < 0.5)),
    "==================="
)
cat(ndr_summary, sep = "\n")
writeLines(ndr_summary, "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu_training/ndr_summary.txt")

# Write output 
#if quant = TRUE use
#writeBambuOutput(se, path = "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu/")
writeBambuOutput(se, path = "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu_training/")
#if quant = FALSE use
#writeToGTF(se, "/nfs/scratch/pierotmi/projects/black_coral/rna/output/bambu/output.gtf")


#############
# For further details on the output see: https://github.com/GoekeLab/bambu/issues/467