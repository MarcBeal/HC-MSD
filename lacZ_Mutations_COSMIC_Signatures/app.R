#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")

library(ggplot2)
library(reshape)
library(shiny)
library(deconstructSigs)

lacZ <- "ATGACCATGATTACGGATTCACTGGAATTCCCGGGGATCCCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAACTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCATCATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTACACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATTTGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGATGAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGCCAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAAGCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACCAGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGATGGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTGATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGCGTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAGTGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCATCTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAACCGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACCCGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAAGCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCTCACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGATGGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCGGCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGATTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGATCTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGCGGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATCAGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAAGAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGGAGCCCGTCAGTATCGGCGGAATTACAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAAAAATAATAATAACCGGGCAGGCCATGTCTGCCCGTATTTCGCGTAAGGAAATCCATTATGT"

################################################################################
################################################################################
##########################Prepare Data, Compare to Sigs#########################
################################################################################
################################################################################

compareToSigs <- function(sigInput, sigVersion=3) {

df_vr <- list()

errorInRef <- vector()

for (i in 1:length(levels(sigInput$Group))) {

	currentGroup <- sigInput[sigInput$Group==levels(sigInput$Group)[i],]
	
	#code to check for discrepancies between Ref and Actual Ref
	#for (i in 1:nrow(currentGroup)) { if (currentGroup$Ref[i] != substr(lacZ, currentGroup$Pos[i], currentGroup$Pos[i])) {print(currentGroup$Pos[i])} }

	df_vr[[i]] <- matrix(nrow=1, ncol=ncol(signatures.cosmic))
	df_vr[[i]][1,] <- rep(0, 96)
	colnames(df_vr[[i]]) <- colnames(signatures.cosmic)
	rownames(df_vr[[i]]) <- levels(sigInput$Group)[i]

	for (j in 1:nrow(currentGroup)){
	
		if(as.character(currentGroup$Ref[j]) != substr(lacZ, currentGroup$Pos[j], currentGroup$Pos[j])) { errorInRef <- c(errorInRef, paste(levels(sigInput$Group)[i], ", Position:", currentGroup$Pos[j], ", Ref:", currentGroup$Ref[j], sep="")) }
	
		if (as.character(currentGroup$Ref[j]) == substr(lacZ, currentGroup$Pos[j], currentGroup$Pos[j]) & currentGroup$Pos[j] != 1) {
			
			if (currentGroup$Ref[j] == "C") {
				if (currentGroup$Alt[j] == "A") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACA") {
					df_vr[[i]][1] = df_vr[[i]][1] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACC") {
					df_vr[[i]][2] = df_vr[[i]][2] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACG") {
					df_vr[[i]][3] = df_vr[[i]][3] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACT") {
					df_vr[[i]][4] = df_vr[[i]][4] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCA") {
					df_vr[[i]][5] = df_vr[[i]][5] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCC") {
					df_vr[[i]][6] = df_vr[[i]][6] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCG") {
					df_vr[[i]][7] = df_vr[[i]][7] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCT") {
					df_vr[[i]][8] = df_vr[[i]][8] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCA") {
					df_vr[[i]][9] = df_vr[[i]][9] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCC") {
					df_vr[[i]][10] = df_vr[[i]][10] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCG") {
					df_vr[[i]][11] = df_vr[[i]][11] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCT") {
					df_vr[[i]][12] = df_vr[[i]][12] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCA") {
					df_vr[[i]][13] = df_vr[[i]][13] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCC") {
					df_vr[[i]][14] = df_vr[[i]][14] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCG") {
					df_vr[[i]][15] = df_vr[[i]][15] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCT") {
					df_vr[[i]][16] = df_vr[[i]][16] + 1
					}					
				}
				if (currentGroup$Alt[j] == "G") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACA") {
					df_vr[[i]][17] = df_vr[[i]][17] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACC") {
					df_vr[[i]][18] = df_vr[[i]][18] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACG") {
					df_vr[[i]][19] = df_vr[[i]][19] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACT") {
					df_vr[[i]][20] = df_vr[[i]][20] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCA") {
					df_vr[[i]][21] = df_vr[[i]][21] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCC") {
					df_vr[[i]][22] = df_vr[[i]][22] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCG") {
					df_vr[[i]][23] = df_vr[[i]][23] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCT") {
					df_vr[[i]][24] = df_vr[[i]][24] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCA") {
					df_vr[[i]][25] = df_vr[[i]][25] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCC") {
					df_vr[[i]][26] = df_vr[[i]][26] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCG") {
					df_vr[[i]][27] = df_vr[[i]][27] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCT") {
					df_vr[[i]][28] = df_vr[[i]][28] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCA") {
					df_vr[[i]][29] = df_vr[[i]][29] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCC") {
					df_vr[[i]][30] = df_vr[[i]][30] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCG") {
					df_vr[[i]][31] = df_vr[[i]][31] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCT") {
					df_vr[[i]][32] = df_vr[[i]][32] + 1
					}					
				}
				if (currentGroup$Alt[j] == "T") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACA") {
					df_vr[[i]][33] = df_vr[[i]][33] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACC") {
					df_vr[[i]][34] = df_vr[[i]][34] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACG") {
					df_vr[[i]][35] = df_vr[[i]][35] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ACT") {
					df_vr[[i]][36] = df_vr[[i]][36] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCA") {
					df_vr[[i]][37] = df_vr[[i]][37] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCC") {
					df_vr[[i]][38] = df_vr[[i]][38] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCG") {
					df_vr[[i]][39] = df_vr[[i]][39] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CCT") {
					df_vr[[i]][40] = df_vr[[i]][40] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCA") {
					df_vr[[i]][41] = df_vr[[i]][41] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCC") {
					df_vr[[i]][42] = df_vr[[i]][42] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCG") {
					df_vr[[i]][43] = df_vr[[i]][43] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GCT") {
					df_vr[[i]][44] = df_vr[[i]][44] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCA") {
					df_vr[[i]][45] = df_vr[[i]][45] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCC") {
					df_vr[[i]][46] = df_vr[[i]][46] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCG") {
					df_vr[[i]][47] = df_vr[[i]][47] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TCT") {
					df_vr[[i]][48] = df_vr[[i]][48] + 1
					}					
				}				
			}
			if (currentGroup$Ref[j] == "T") {
				if (currentGroup$Alt[j] == "A") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATA") {
					df_vr[[i]][49] = df_vr[[i]][49] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATC") {
					df_vr[[i]][50] = df_vr[[i]][50] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATG") {
					df_vr[[i]][51] = df_vr[[i]][51] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATT") {
					df_vr[[i]][52] = df_vr[[i]][52] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTA") {
					df_vr[[i]][53] = df_vr[[i]][53] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTC") {
					df_vr[[i]][54] = df_vr[[i]][54] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTG") {
					df_vr[[i]][55] = df_vr[[i]][55] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTT") {
					df_vr[[i]][56] = df_vr[[i]][56] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTA") {
					df_vr[[i]][57] = df_vr[[i]][57] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTC") {
					df_vr[[i]][58] = df_vr[[i]][58] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTG") {
					df_vr[[i]][59] = df_vr[[i]][59] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTT") {
					df_vr[[i]][60] = df_vr[[i]][60] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTA") {
					df_vr[[i]][61] = df_vr[[i]][61] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTC") {
					df_vr[[i]][62] = df_vr[[i]][62] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTG") {
					df_vr[[i]][63] = df_vr[[i]][63] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTT") {
					df_vr[[i]][64] = df_vr[[i]][64] + 1
					}					
				}
				if (currentGroup$Alt[j] == "C") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATA") {
					df_vr[[i]][65] = df_vr[[i]][65] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATC") {
					df_vr[[i]][66] = df_vr[[i]][66] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATG") {
					df_vr[[i]][67] = df_vr[[i]][67] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATT") {
					df_vr[[i]][68] = df_vr[[i]][68] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTA") {
					df_vr[[i]][69] = df_vr[[i]][69] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTC") {
					df_vr[[i]][70] = df_vr[[i]][70] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTG") {
					df_vr[[i]][71] = df_vr[[i]][71] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTT") {
					df_vr[[i]][72] = df_vr[[i]][72] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTA") {
					df_vr[[i]][73] = df_vr[[i]][73] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTC") {
					df_vr[[i]][74] = df_vr[[i]][74] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTG") {
					df_vr[[i]][75] = df_vr[[i]][75] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTT") {
					df_vr[[i]][76] = df_vr[[i]][76] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTA") {
					df_vr[[i]][77] = df_vr[[i]][77] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTC") {
					df_vr[[i]][78] = df_vr[[i]][78] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTG") {
					df_vr[[i]][79] = df_vr[[i]][79] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTT") {
					df_vr[[i]][80] = df_vr[[i]][80] + 1
					}					
				}
				if (currentGroup$Alt[j] == "G") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATA") {
					df_vr[[i]][81] = df_vr[[i]][81] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATC") {
					df_vr[[i]][82] = df_vr[[i]][82] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATG") {
					df_vr[[i]][83] = df_vr[[i]][83] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "ATT") {
					df_vr[[i]][84] = df_vr[[i]][84] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTA") {
					df_vr[[i]][85] = df_vr[[i]][85] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTC") {
					df_vr[[i]][86] = df_vr[[i]][86] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTG") {
					df_vr[[i]][87] = df_vr[[i]][87] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CTT") {
					df_vr[[i]][88] = df_vr[[i]][88] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTA") {
					df_vr[[i]][89] = df_vr[[i]][89] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTC") {
					df_vr[[i]][90] = df_vr[[i]][90] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTG") {
					df_vr[[i]][91] = df_vr[[i]][91] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GTT") {
					df_vr[[i]][92] = df_vr[[i]][92] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTA") {
					df_vr[[i]][93] = df_vr[[i]][93] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTC") {
					df_vr[[i]][94] = df_vr[[i]][94] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTG") {
					df_vr[[i]][95] = df_vr[[i]][95] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TTT") {
					df_vr[[i]][96] = df_vr[[i]][96] + 1
					}					
				}
			
			}
			if (currentGroup$Ref[j] == "G") {
				if (currentGroup$Alt[j] == "T") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGT") {
					df_vr[[i]][1] = df_vr[[i]][1] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGT") {
					df_vr[[i]][2] = df_vr[[i]][2] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGT") {
					df_vr[[i]][3] = df_vr[[i]][3] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGT") {
					df_vr[[i]][4] = df_vr[[i]][4] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGG") {
					df_vr[[i]][5] = df_vr[[i]][5] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGG") {
					df_vr[[i]][6] = df_vr[[i]][6] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGG") {
					df_vr[[i]][7] = df_vr[[i]][7] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGG") {
					df_vr[[i]][8] = df_vr[[i]][8] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGC") {
					df_vr[[i]][9] = df_vr[[i]][9] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGC") {
					df_vr[[i]][10] = df_vr[[i]][10] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGC") {
					df_vr[[i]][11] = df_vr[[i]][11] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGC") {
					df_vr[[i]][12] = df_vr[[i]][12] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGA") {
					df_vr[[i]][13] = df_vr[[i]][13] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGA") {
					df_vr[[i]][14] = df_vr[[i]][14] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGA") {
					df_vr[[i]][15] = df_vr[[i]][15] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGA") {
					df_vr[[i]][16] = df_vr[[i]][16] + 1
					}					
				}
				if (currentGroup$Alt[j] == "C") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGT") {
					df_vr[[i]][17] = df_vr[[i]][17] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGT") {
					df_vr[[i]][18] = df_vr[[i]][18] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGT") {
					df_vr[[i]][19] = df_vr[[i]][19] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGT") {
					df_vr[[i]][20] = df_vr[[i]][20] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGG") {
					df_vr[[i]][21] = df_vr[[i]][21] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGG") {
					df_vr[[i]][22] = df_vr[[i]][22] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGG") {
					df_vr[[i]][23] = df_vr[[i]][23] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGG") {
					df_vr[[i]][24] = df_vr[[i]][24] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGC") {
					df_vr[[i]][25] = df_vr[[i]][25] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGC") {
					df_vr[[i]][26] = df_vr[[i]][26] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGC") {
					df_vr[[i]][27] = df_vr[[i]][27] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGC") {
					df_vr[[i]][28] = df_vr[[i]][28] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGA") {
					df_vr[[i]][29] = df_vr[[i]][29] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGA") {
					df_vr[[i]][30] = df_vr[[i]][30] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGA") {
					df_vr[[i]][31] = df_vr[[i]][31] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGA") {
					df_vr[[i]][32] = df_vr[[i]][32] + 1
					}					
				}
				if (currentGroup$Alt[j] == "A") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGT") {
					df_vr[[i]][33] = df_vr[[i]][33] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGT") {
					df_vr[[i]][34] = df_vr[[i]][34] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGT") {
					df_vr[[i]][35] = df_vr[[i]][35] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGT") {
					df_vr[[i]][36] = df_vr[[i]][36] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGG") {
					df_vr[[i]][37] = df_vr[[i]][37] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGG") {
					df_vr[[i]][38] = df_vr[[i]][38] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGG") {
					df_vr[[i]][39] = df_vr[[i]][39] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGG") {
					df_vr[[i]][40] = df_vr[[i]][40] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGC") {
					df_vr[[i]][41] = df_vr[[i]][41] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGC") {
					df_vr[[i]][42] = df_vr[[i]][42] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGC") {
					df_vr[[i]][43] = df_vr[[i]][43] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGC") {
					df_vr[[i]][44] = df_vr[[i]][44] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TGA") {
					df_vr[[i]][45] = df_vr[[i]][45] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GGA") {
					df_vr[[i]][46] = df_vr[[i]][46] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CGA") {
					df_vr[[i]][47] = df_vr[[i]][47] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AGA") {
					df_vr[[i]][48] = df_vr[[i]][48] + 1
					}					
				}				
			}
			if (currentGroup$Ref[j] == "A") {
				if (currentGroup$Alt[j] == "T") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAT") {
					df_vr[[i]][49] = df_vr[[i]][49] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAT") {
					df_vr[[i]][50] = df_vr[[i]][50] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAT") {
					df_vr[[i]][51] = df_vr[[i]][51] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAT") {
					df_vr[[i]][52] = df_vr[[i]][52] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAG") {
					df_vr[[i]][53] = df_vr[[i]][53] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAG") {
					df_vr[[i]][54] = df_vr[[i]][54] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAG") {
					df_vr[[i]][55] = df_vr[[i]][55] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAG") {
					df_vr[[i]][56] = df_vr[[i]][56] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAC") {
					df_vr[[i]][57] = df_vr[[i]][57] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAC") {
					df_vr[[i]][58] = df_vr[[i]][58] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAC") {
					df_vr[[i]][59] = df_vr[[i]][59] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAC") {
					df_vr[[i]][60] = df_vr[[i]][60] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAA") {
					df_vr[[i]][61] = df_vr[[i]][61] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAA") {
					df_vr[[i]][62] = df_vr[[i]][62] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAA") {
					df_vr[[i]][63] = df_vr[[i]][63] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAA") {
					df_vr[[i]][64] = df_vr[[i]][64] + 1
					}					
				}
				if (currentGroup$Alt[j] == "G") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAT") {
					df_vr[[i]][65] = df_vr[[i]][65] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAT") {
					df_vr[[i]][66] = df_vr[[i]][66] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAT") {
					df_vr[[i]][67] = df_vr[[i]][67] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAT") {
					df_vr[[i]][68] = df_vr[[i]][68] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAG") {
					df_vr[[i]][69] = df_vr[[i]][69] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAG") {
					df_vr[[i]][70] = df_vr[[i]][70] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAG") {
					df_vr[[i]][71] = df_vr[[i]][71] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAG") {
					df_vr[[i]][72] = df_vr[[i]][72] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAC") {
					df_vr[[i]][73] = df_vr[[i]][73] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAC") {
					df_vr[[i]][74] = df_vr[[i]][74] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAC") {
					df_vr[[i]][75] = df_vr[[i]][75] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAC") {
					df_vr[[i]][76] = df_vr[[i]][76] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAA") {
					df_vr[[i]][77] = df_vr[[i]][77] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAA") {
					df_vr[[i]][78] = df_vr[[i]][78] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAA") {
					df_vr[[i]][79] = df_vr[[i]][79] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAA") {
					df_vr[[i]][80] = df_vr[[i]][80] + 1
					}					
				}
				if (currentGroup$Alt[j] == "C") {
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAT") {
					df_vr[[i]][81] = df_vr[[i]][81] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAT") {
					df_vr[[i]][82] = df_vr[[i]][82] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAT") {
					df_vr[[i]][83] = df_vr[[i]][83] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAT") {
					df_vr[[i]][84] = df_vr[[i]][84] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAG") {
					df_vr[[i]][85] = df_vr[[i]][85] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAG") {
					df_vr[[i]][86] = df_vr[[i]][86] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAG") {
					df_vr[[i]][87] = df_vr[[i]][87] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAG") {
					df_vr[[i]][88] = df_vr[[i]][88] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAC") {
					df_vr[[i]][89] = df_vr[[i]][89] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAC") {
					df_vr[[i]][90] = df_vr[[i]][90] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAC") {
					df_vr[[i]][91] = df_vr[[i]][91] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAC") {
					df_vr[[i]][92] = df_vr[[i]][92] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "TAA") {
					df_vr[[i]][93] = df_vr[[i]][93] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "GAA") {
					df_vr[[i]][94] = df_vr[[i]][94] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "CAA") {
					df_vr[[i]][95] = df_vr[[i]][95] + 1
					}
					if (substr(lacZ, currentGroup$Pos[j] - 1, currentGroup$Pos[j] + 1) == "AAA") {
					df_vr[[i]][96] = df_vr[[i]][96] + 1
					}					
				}
			
			}
			
		}	
	}

	df_vr[[i]] <- df_vr[[i]] / sum(df_vr[[i]])
	df_vr[[i]] <- data.frame(df_vr[[i]])
	
	

}

#Data ready for whichSignatures

	for (i in 1:length(levels(sigInput$Group))) {
		colnames(df_vr[[i]]) <- colnames(signatures.cosmic)
	}

	if (sigVersion == 3) {
		lacZ_cosmic_signatures <- as.data.frame(t(read.table("data/lacZ_corrected_COSMIC_signatures_v3_with_control.txt", sep="\t", header=TRUE)))
	}
	if (sigVersion == 2) {
		lacZ_cosmic_signatures <- as.data.frame(t(read.table("data/lacZ_corrected_COSMIC_signatures_v2_with_control.txt", sep="\t", header=TRUE)))
	}

	colnames(lacZ_cosmic_signatures) <- colnames(signatures.cosmic)


	plotSigs <- list()

	for (i in 1:length(levels(sigInput$Group))) {
		plotSigs[[i]] <- whichSignatures(tumor.ref=df_vr[[i]], signatures.ref=lacZ_cosmic_signatures, sample.id=rownames(df_vr[[i]]))
	}

	#Extract the signature data for each study group

	allSigs <- matrix(ncol = nrow(lacZ_cosmic_signatures), nrow = length(levels(sigInput$Group)))

	vcfNames <- vector()

	for (i in 1:length(levels(sigInput$Group))) {
		allSigs[i,] <- as.numeric(plotSigs[[i]]$weights)
		vcfNames <- c(vcfNames, rownames(plotSigs[[i]]$weights))
	}

	rownames(allSigs) <- vcfNames
	colnames(allSigs) <- names(plotSigs[[1]]$weights)


	resultsToReturn <- list()

	resultsToReturn$signatureSummary <- allSigs
	resultsToReturn$plotSigs <- plotSigs
	resultsToReturn$lacZ_cosmic_signatures <- lacZ_cosmic_signatures
	resultsToReturn$errorInRef <- errorInRef 

	return(resultsToReturn)

}



################################################################################
################################################################################
############################Plotting Function Begins############################
################################################################################
################################################################################


#inData comes from the data produced from the previous function
#vcfNumbers are the corresponding sample groups that you want to visualize that can be plotted when lacZSigs is False
#If you would like to plot existing signatures, set lacZSigs to TRUE and give the signature numbers as a vector (example: c(3, 4, 18))


plotSignatureResults <- function(inData, vcfNumbers, lacZSigs=F, whichSigs) {


if (lacZSigs == T) {
	lacZ_cosmic_signatures <- inData$lacZ_cosmic_signatures
	lacZ_cosmic_signatures$Signature <- row.names(lacZ_cosmic_signatures)
	lacZ_cosmic_signatures_melted <- melt(lacZ_cosmic_signatures, id.vars = "Signature" )
	tempmuts <-gsub("[A-Z]\\[", "", lacZ_cosmic_signatures_melted$variable)
	lacZ_cosmic_signatures_melted$muttype <- gsub("\\][A-Z]", "", tempmuts)
	lacZ_cosmic_signatures_melted$Signature_factor <- factor(lacZ_cosmic_signatures_melted$Signature, c(unique(lacZ_cosmic_signatures_melted$Signature)))

	listSigs <- vector()
	for (i in 1:length(whichSigs)) {
		listSigs <- c(listSigs, paste("Signature.", whichSigs[i], sep=""))
	}

	ggplot(lacZ_cosmic_signatures_melted[lacZ_cosmic_signatures_melted$Signature_factor %in% listSigs,], aes(x=variable, y=value)) + geom_bar(stat="identity", aes(fill=muttype), width=0.8) + scale_fill_manual(values = c("#1ebff0", "#000000", "#e62725", "#cbcacb", "#a1cf64", "#edc8c5")) + facet_grid(Signature_factor ~ muttype, scales="free_x") + theme(panel.spacing.x=unit(0, "lines"), panel.grid.major = element_line(colour="grey", size=0.5), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks = NULL) + scale_x_discrete(breaks=NULL) + coord_cartesian(ylim=c(0,0.4))

}

if (lacZSigs == F) {
	vcfsToPlot <- matrix(ncol=96)
	colnames(vcfsToPlot) <- colnames(inData$plotSigs[[1]]$tumor)
	
	for (i in 1:length(vcfNumbers)) {	
		vcfsToPlot <- rbind(vcfsToPlot, inData$plotSigs[[vcfNumbers[i]]]$tumor)
		}
	vcfsToPlot <- as.data.frame(vcfsToPlot[-1,])

	vcfsToPlot$tissue <- rownames(vcfsToPlot)
	vcfsToPlot_melted <- melt(vcfsToPlot, id.vars = "tissue" )
	tempmuts <-gsub("[A-Z]\\[", "", vcfsToPlot_melted$variable)
	vcfsToPlot_melted$muttype <- gsub("\\][A-Z]", "", tempmuts)
	vcfsToPlot_melted$Tissue_factor <- factor(vcfsToPlot_melted$tissue, c(unique(vcfsToPlot_melted$tissue)))

	

	listvcfs <- vector()
	for (i in 1:length(vcfNumbers)) {
		listvcfs <- c(listvcfs, rownames(vcfsToPlot)[i])		
	}

	

	if (length(vcfNumbers) < 2) { print("vcfNumbers needs to include more than 1 data set to work properly") } else {
	ggplot(vcfsToPlot_melted[vcfsToPlot_melted$Tissue_factor %in% listvcfs,], aes(x=variable, y=value)) + geom_bar(stat="identity", aes(fill=muttype), width=0.8) + scale_fill_manual(values = c("#1ebff0", "#000000", "#e62725", "#cbcacb", "#a1cf64", "#edc8c5")) + facet_grid(Tissue_factor ~ muttype, scales="free_x") + theme(panel.spacing.x=unit(0, "lines"), panel.grid.major = element_line(colour="grey", size=0.5), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks = NULL) + scale_x_discrete(breaks=NULL) + coord_cartesian(ylim=c(0,0.4))	
	}
}



}

################################################################################
################################################################################
###################Reconstruct & Similarity Comparison Starts###################
################################################################################
################################################################################

reconstructSig <- function(inData, simMetric) {

reconstructed <- list()
p.value <- vector()
correlation <- vector()

	#Calculate Pearson Correlation and Cosine Similarity
	for (i in 1:nrow(inData$signatureSummary)) {
		reconstructed[[i]] <- inData$lacZ_cosmic_signatures * inData$signatureSummary[i,]
		
		#Below calculates Pearson Correlation
		result <- cor.test(as.vector(colSums(reconstructed[[i]])), as.vector(inData$plotSigs[[i]]$tumor), "two.sided", "pearson")
		p.value <- c(p.value, result$p.value)
		
		if (simMetric == "Pearson Correlation") {
				correlation <- c(correlation, result$estimate)
		}
		
		#Below reports Cosine Similarity
		mutPattern1 <- as.vector(colSums(reconstructed[[i]]))
		mutPattern2 <- as.vector(inData$plotSigs[[i]]$tumor)
		
		if (simMetric == "Cosine Similarity") {
			correlation <- c(correlation, sum(mutPattern1*mutPattern2)/sqrt(sum(mutPattern1^2)*sum(mutPattern2^2)))
		}
	}

resultsSummary <- rbind(p.value, correlation) #note: p-value corresponds to Pearson correlation and is not reported in output
colnames(resultsSummary) <- rownames(inData$signatureSummary)

return(resultsSummary)

}

################################################################################
################################################################################
############Correlation Between Ind. Sigs. and Mutation Data. Starts############
################################################################################
################################################################################

corTestSig <- function(inData, simMetric) {

group <- list()

	for (i in 1:nrow(inData$signatureSummary)) {

	group[[i]] <- list()
	group[[i]]$groupName <- row.names(inData$signatureSummary)[i]
	
	tested <- vector()
	#p.value <- vector()
	correlation <- vector()
		
	corSigs <- as.integer(which(inData$signatureSummary[i,]>0))

		for (j in 1:length(corSigs)) {
			#Below records Pearson Correlation
			result <- cor.test(as.numeric(inData$lacZ_cosmic_signatures[corSigs[j],]), as.numeric(inData$plotSigs[[i]]$tumor), "two.sided", "pearson")		
			tested <- c(tested, colnames(inData$signatureSummary)[corSigs][j])
						
			if (simMetric == "Pearson Correlation") {
				correlation <- c(correlation, result$estimate)
			} else if (simMetric == "Cosine Similarity") {
				#Below reports Cosine Similarity - Updated analysis
				mutPattern1 <- as.numeric(inData$lacZ_cosmic_signatures[corSigs[j],])
				mutPattern2 <- as.numeric(inData$plotSigs[[i]]$tumor)
				correlation <- c(correlation, sum(mutPattern1*mutPattern2)/sqrt(sum(mutPattern1^2)*sum(mutPattern2^2)))
			}
		
		}

	group[[i]]$result <- rbind(tested, correlation)

	}

return(group)

}

################################################################################
################################################################################
#########################Generate Random Mutations##############################
################################################################################
################################################################################


randomMutationGenerator <- function(numMutations=1000, groupName="Random") {

Group <- rep(groupName, numMutations)
Position <- sample(2:3096, numMutations, replace=T)
Ref <- vector()
Alt <- vector()




for(i in 1:length(Position)) {
	
	Ref <- c(Ref, substr(lacZ, Position[i], Position[i]))
	
	mutation <- sample(c("A", "C", "G", "T"), 1)
	
	while (mutation == Ref[i]) {
		mutation <- sample(c("A", "C", "G", "T"), 1)
	}
	
	Alt <- c(Alt, mutation)
	
}

return(data.frame(cbind(Group, Position, Ref, Alt)))

}

################################################################################
################################################################################
#################################GUI Begins#####################################
################################################################################
################################################################################



ui <- fluidPage(

  titlePanel(" ", "Transgenic-COSMIC Analysis Tool"),
  
	sidebarPanel(

  tags$table(
    tags$tr(
      #tags$td(img(src="beta-gal.png", height=80)),
      tags$td(h3("Compare ", em("lacZ"), " Mutations With Cosmic Signatures"))
    )
  ),
	
  fileInput("file1", "Choose Tab-Delimited File",
                multiple = FALSE,
                accept = c("text/tab-separated-values")),

	helpText("Please upload file with exact columns: Group, Position, Ref, Alt. After upload is completed click on Analyze and please wait a moment for analysis to complete. Do not hit Analyze again or any of the Download buttons until analysis is completed. This will be evident by the appearance of a summary message on the right. Downloads may take a moment to begin. Currently, a minimum of two groups are required to run program. Graph Y-axis maxes out at 0.4"),

	radioButtons("cosmicVersion", "Choose COSMIC Mutational Signatures Version", c("Version 3", "Version 2"), selected="Version 3"),
	
	radioButtons("simMetric", "Choose Similarity Metric for Data Filtering", c("Cosine Similarity", "Pearson Correlation"), selected="Cosine Similarity"),

	actionButton("analyze", "Analyze"),

	hr(),

	downloadButton("downloadRawData", "Download Raw Summary"),

	downloadButton("downloadFinalData", "Download Final Summary"),

	downloadButton("downloadCorResults", "Download Similarity Metrics"),

	downloadButton("downloadGraphs", "Download Graphs"),

	downloadButton("downloadGraphData", "Download Graph Data"),

	hr(),

	actionLink("linkLacZ", "View Reference FASTA"),
	
	verbatimTextOutput("displayLacZ"),
	
	hr(),
	
	tags$h3("Random Mutation Generator Tool"),
	
	numericInput("numMuts", "Enter Number of Mutations Desired", 1000),
	
	numericInput("numRandGroups", "Enter Number of Random Groups", 2),
	
	helpText("Table created can be read by comparison tool above."),
	
	downloadButton("downloadRand", "Create and Download Random Mutations")

	),

	mainPanel(

	verbatimTextOutput("summaryText"),

	plotOutput("signaturePlots", height=2000)

	)

)

server <- function(input, output) {

	analyzeData <- eventReactive(input$analyze, {

	withProgress(message = "Analyzing Raw Data", value = 0, {	

		inFile = input$file1
		sigInput <- read.table(inFile$datapath, header = T, sep="\t")
		
		incProgress(0.5)
		
		if (input$cosmicVersion == "Version 3") { rawResults <- compareToSigs(sigInput, 3) }
		if (input$cosmicVersion == "Version 2") { rawResults <- compareToSigs(sigInput, 2) }

		incProgress(0.5)
		Sys.sleep(0.1)

		rawResults

		})
	})

	catLacZ <- eventReactive(input$linkLacZ, { 
		cat(">lacZ reference\n")
		cat(lacZ)
	})

	output$displayLacZ <- renderPrint({ catLacZ() })

	output$downloadRand <- downloadHandler(
	filename = function() {
	"RandomMutations.txt"
	},
	content = function(file) {
	
	tempTable <- randomMutationGenerator(input$numMuts, "Random1")
	outputCols <- colnames(tempTable)
	
	for (i in 2:input$numRandGroups) {
		tempTable <- rbind(tempTable, randomMutationGenerator(input$numMuts, paste("Random", i, sep="")))
	}
	
	colnames(tempTable) <- outputCols
	
	write.table(tempTable, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	})

	output$signaturePlots <- renderPlot({
		plotSignatureResults(analyzeData(), 1:nrow(analyzeData()$signatureSummary)) 
	})

	###Remove zero sigs, merge summary with pearson/cosine (from reconstructSigs), print out raw summary
	rawSummaryOutput <- reactive({

		rawSummary <- analyzeData()$signatureSummary[, colSums(analyzeData()$signatureSummary != 0) > 0]

		Residual <- 1 - t(t(rowSums(rawSummary)))

		rawSummary <- cbind(rawSummary, Residual)
	
		colnames(rawSummary)[ncol(rawSummary)] = "Residual"

		rawSummary <- cbind(rawSummary, t(t(reconstructSig(analyzeData(), input$simMetric)[2,])))

		if (input$simMetric == "Pearson Correlation") {
			colnames(rawSummary)[ncol(rawSummary)] = "reconstructSigPearsonCoefficient"
		} else if (input$simMetric == "Cosine Similarity") {
			colnames(rawSummary)[ncol(rawSummary)] = "reconstructSigCosineSimilarity"
		}
		
		rawSummary <- t(rawSummary)
		
	})

	output$downloadRawData <- downloadHandler(
    	filename = function() {
      	"Signature_raw_summary.txt"
    	},
    	content = function(file) {
      	write.table(rawSummaryOutput(), file, sep="\t", col.names=NA)
    	})


	###Write all correlations

	output$downloadCorResults <- downloadHandler(
    	filename = function() {
      	"Signature_similarity_between_sigs_and_mut_profiles.txt"
		},
    	content = function(file) {
      	correlationResults <- corTestSig(analyzeData(), input$simMetric)
	sink(file)
    	print(correlationResults)
	sink()
	})


	output$downloadFinalData <- downloadHandler(
    	filename = function() {
      	"Signature_final_summary_above_thresholds.txt"
    	},
    	content = function(file) {

			#Find largest residual and filter out signatures with contributions below	
			rawSummary <- analyzeData()$signatureSummary[, colSums(analyzeData()$signatureSummary != 0) > 0]
			Residual <- 1 - t(t(rowSums(rawSummary)))      	
			largestResidual <- max(Residual)

			finalSummary <- analyzeData()$signatureSummary[, colSums(analyzeData()$signatureSummary != 0) > 0]
			finalSummary <- ifelse(finalSummary<largestResidual,0,finalSummary)

			#now remove anything with cosine similarity or pearson correlation below 0.5
			correlationResults <- corTestSig(analyzeData(), input$simMetric)
			for(i in 1:nrow(finalSummary)) {
				tempData <- t(correlationResults[[i]]$result)
				nonImportantSigs <- as.vector(t(tempData[tempData[,"correlation"]<=0.5,1]))
				finalSummary[i,][nonImportantSigs] <- 0
			}


			#write final results that pass the thresholds
			finalSummary <- cbind(finalSummary, t(t(reconstructSig(analyzeData(), input$simMetric)["correlation",])))
			
			if (input$simMetric == "Pearson Correlation") {
				colnames(finalSummary)[ncol(finalSummary)] = "reconstructSigPearsonCoefficient"
			} else if (input$simMetric == "Cosine Similarity") {
				colnames(finalSummary)[ncol(finalSummary)] = "reconstructSigCosineSimilarity"
			}
						
			finalSummary <- t(finalSummary)

			write.table(finalSummary, file, sep="\t", col.names=NA)
	})

	output$downloadGraphs <- downloadHandler(
    	filename = function() {
      	"Mutational_Profiles.pdf"
    	},
    	content = function(file) {
	pdf(file)
	print(plotSignatureResults(analyzeData(), 1:nrow(analyzeData()$signatureSummary)))
	dev.off()
 	})

	output$downloadGraphData <- downloadHandler(
	filename = function() {
	"Mutational_Profiles_Data.txt"
	},
	content = function(file) {

	graphData <- analyzeData()$plotSig
	
	tableGraphData <- matrix(nrow=length(graphData), ncol=ncol(graphData[[1]]$tumor))

	rownames(tableGraphData) <- 1:length(graphData)
	
	for (i in 1:length(graphData)) {
		tableGraphData[i,] <- graphData[[i]]$tumor
		rownames(tableGraphData)[i] <- rownames(graphData[[i]]$tumor)[1]
	}

	colnames(tableGraphData) <- colnames(graphData[[1]]$tumor)
	
	write.table(as.data.frame(tableGraphData), file, sep="\t", col.names=NA)
	})


	output$summaryText <- renderPrint({
	rawSummary <- analyzeData()$signatureSummary[, colSums(analyzeData()$signatureSummary != 0) > 0]
	Residual <- 1 - t(t(rowSums(rawSummary)))      	
	largestResidual <- max(Residual)
	cat(paste("Analysis was completed on ", Sys.time(), ".\n", sep=""))
	cat("\n")
	
	if (length(analyzeData()$errorInRef) > 0) {
	cat("Warning! The following mutations were excluded from analysis\nbecause the reference nucleotide reported did not match the\nreference sequence in lacZ:\n\t")
	cat(paste(analyzeData()$errorInRef, collapse="\n\t"))
	cat("\nMake sure positions match lacZ reference.\n\n")
	}

	cat("The raw results from signature deconstruction can\n")
	cat("be found in Signature_raw_summary.txt by clicking\n")
	cat("Download Raw Summary.\n")
	cat("\n")
	cat("The final results can be found in\n")
	cat("Signature_final_summary_above_thresholds.txt\n")
	cat("by clicking Download Final Summary.\n")
	cat("\n")
	cat("Signatures that had a contribution below the\n")
	cat(paste("largest residual of ", largestResidual, " were\n", sep=""))
	cat("removed and reported as 0.\n")
	cat("\n")
	if (input$simMetric == "Pearson Correlation") {
		cat("Furthermore, if the Pearson Correlation between\n")
	} else if (input$simMetric == "Cosine Similarity") {
		cat("Furthermore, if the Cosine Similarity between\n")
	}
	cat("the signature and mutation data was below 0.5,\n")
	cat("that signature was also removed and reported as 0.\n")
	cat("\n")
	cat("The raw P.C. information can be found in\n")
	cat("Signature_similarities_between_sigs_and_mut_profiles.txt\n")
	cat("by clicking Download Similarity Metrics.\n")
	cat("\n")
	cat("Best of luck with your data interpretation!")
	})
	

}

shinyApp(ui, server)
