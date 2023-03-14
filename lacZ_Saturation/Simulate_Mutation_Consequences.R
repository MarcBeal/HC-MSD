###This code will simulate the consequences of all possible mutations in lacZ

library("Biostrings")

lacZ <- "ATGACCATGATTACGGATTCACTGGAATTCCCGGGGATCCCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAACTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCATCATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTACACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATTTGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGATGAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGCCAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAAGCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACCAGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGATGGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTGATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGCGTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAGTGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCATCTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAACCGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACCCGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAAGCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCTCACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGATGGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCGGCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGATTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGATCTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGCGGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATCAGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAAGAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGGAGCCCGTCAGTATCGGCGGAATTACAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAAAAATAATAATAA"

codons <- vector()
startPos <- vector()

#extract the codons and their start position
i = 1
while (i <= nchar(lacZ)) {

codons <- c(codons,substr(lacZ, i, i+2))
startPos <- c(startPos, i)
i = i + 3

}


mutPos <- vector()
mutAlt <- vector()
mutResult <- vector()
refVec <- vector()


#create mutants and see if they are missense, silent, or nonsense
#skip the first codon because all 9 mutations are loss of start codon (already recorded in excel)

j = 2
while (j <= length(codons)) {

	currentNuc <- startPos[j]
	refCodon <- codons[j]
	
	#for each mutation at first position of codon add to the vector the position
	mutPos <- c(mutPos, currentNuc, currentNuc, currentNuc)
	refNuc <- substr(refCodon, 1, 1)

	#first position, make mutations and determine outcome

	mut1 <- codons[j]
	mut2 <- codons[j]
	mut3 <- codons[j]
	

	if (refNuc == "A") {

		mutAlt <- c(mutAlt, "T", "C", "G")
		refVec <- c(refVec, "A", "A", "A")
		
		substr(mut1, 1, 1) = "T"
		substr(mut2, 1, 1) = "C"
		substr(mut3, 1, 1) = "G"
	}

	if (refNuc == "T") {

		mutAlt <- c(mutAlt, "A", "C", "G")
		refVec <- c(refVec, "T", "T", "T")
		
		substr(mut1, 1, 1) = "A"
		substr(mut2, 1, 1) = "C"
		substr(mut3, 1, 1) = "G"
	}

	if (refNuc == "C") {

		mutAlt <- c(mutAlt, "T", "A", "G")
		refVec <- c(refVec, "C", "C", "C")
		
		substr(mut1, 1, 1) = "T"
		substr(mut2, 1, 1) = "A"
		substr(mut3, 1, 1) = "G"
	}

	if (refNuc == "G") {
	
		mutAlt <- c(mutAlt, "T", "C", "A")
		refVec <- c(refVec, "G", "G", "G")
		
		substr(mut1, 1, 1) = "T"
		substr(mut2, 1, 1) = "C"
		substr(mut3, 1, 1) = "A"
	}


	if (GENETIC_CODE[[mut1]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut1]] != "*" && GENETIC_CODE[[mut1]] != GENETIC_CODE[[codons[j]]]) {mutResult <-  c(mutResult, "missense")}
	if (GENETIC_CODE[[mut1]] == GENETIC_CODE[[codons[j]]]) {mutResult <-  c(mutResult, "silent")}

	if (GENETIC_CODE[[mut2]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <-  c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut2]] != "*" && GENETIC_CODE[[mut2]] != GENETIC_CODE[[codons[j]]]) {mutResult <-  c(mutResult, "missense")}
	if (GENETIC_CODE[[mut2]] == GENETIC_CODE[[codons[j]]]) {mutResult <-  c(mutResult, "silent")}
	
	if (GENETIC_CODE[[mut3]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut3]] != "*" && GENETIC_CODE[[mut3]] != GENETIC_CODE[[codons[j]]]) {mutResult <-  c(mutResult, "missense")}
	if (GENETIC_CODE[[mut3]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}

	#for each mutation at 2nd position of codon add to the vector the position
	mutPos <- c(mutPos, currentNuc + 1, currentNuc + 1, currentNuc + 1)
	refNuc <- substr(refCodon, 2, 2)

	#second position, make mutations and determine outcome

	mut1 <- codons[j]
	mut2 <- codons[j]
	mut3 <- codons[j]
	

	if (refNuc == "A") {

		mutAlt <- c(mutAlt, "T", "C", "G")
		refVec <- c(refVec, "A", "A", "A")
		
		substr(mut1, 2, 2) = "T"
		substr(mut2, 2, 2) = "C"
		substr(mut3, 2, 2) = "G"
	}

	if (refNuc == "T") {

		mutAlt <- c(mutAlt, "A", "C", "G")
		refVec <- c(refVec, "T", "T", "T")
		
		substr(mut1, 2, 2) = "A"
		substr(mut2, 2, 2) = "C"
		substr(mut3, 2, 2) = "G"
	}

	if (refNuc == "C") {

		mutAlt <- c(mutAlt, "T", "A", "G")
		refVec <- c(refVec, "C", "C", "C")
		
		substr(mut1, 2, 2) = "T"
		substr(mut2, 2, 2) = "A"
		substr(mut3, 2, 2) = "G"
	}

	if (refNuc == "G") {
	
		mutAlt <- c(mutAlt, "T", "C", "A")
		refVec <- c(refVec, "G", "G", "G")
		
		substr(mut1, 2, 2) = "T"
		substr(mut2, 2, 2) = "C"
		substr(mut3, 2, 2) = "A"
	}


	if (GENETIC_CODE[[mut1]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut1]] != "*" && GENETIC_CODE[[mut1]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut1]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}

	if (GENETIC_CODE[[mut2]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut2]] != "*" && GENETIC_CODE[[mut2]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut2]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}
	
	if (GENETIC_CODE[[mut3]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut3]] != "*" && GENETIC_CODE[[mut3]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut3]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}


	#for each mutation at 3rd position of codon add to the vector the position
	mutPos <- c(mutPos, currentNuc + 2, currentNuc + 2, currentNuc + 2)
	refNuc <- substr(refCodon, 3, 3)

	#third position, make mutations and determine outcome

	mut1 <- codons[j]
	mut2 <- codons[j]
	mut3 <- codons[j]
	

	if (refNuc == "A") {

		mutAlt <- c(mutAlt, "T", "C", "G")
		refVec <- c(refVec, "A", "A", "A")
		
		substr(mut1, 3, 3) = "T"
		substr(mut2, 3, 3) = "C"
		substr(mut3, 3, 3) = "G"
	}

	if (refNuc == "T") {

		mutAlt <- c(mutAlt, "A", "C", "G")
		refVec <- c(refVec, "T", "T", "T")
		
		substr(mut1, 3, 3) = "A"
		substr(mut2, 3, 3) = "C"
		substr(mut3, 3, 3) = "G"
	}

	if (refNuc == "C") {

		mutAlt <- c(mutAlt, "T", "A", "G")
		refVec <- c(refVec, "C", "C", "C")
		
		substr(mut1, 3, 3) = "T"
		substr(mut2, 3, 3) = "A"
		substr(mut3, 3, 3) = "G"
	}

	if (refNuc == "G") {
	
		mutAlt <- c(mutAlt, "T", "C", "A")
		refVec <- c(refVec, "G", "G", "G")		

		substr(mut1, 3, 3) = "T"
		substr(mut2, 3, 3) = "C"
		substr(mut3, 3, 3) = "A"
	}


	if (GENETIC_CODE[[mut1]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut1]] != "*" && GENETIC_CODE[[mut1]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut1]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}

	if (GENETIC_CODE[[mut2]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut2]] != "*" && GENETIC_CODE[[mut2]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut2]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}
	
	if (GENETIC_CODE[[mut3]] == "*" && GENETIC_CODE[[codons[j]]] != "*") {mutResult <- c(mutResult, "nonsense")}
	if (GENETIC_CODE[[mut3]] != "*" && GENETIC_CODE[[mut3]] != GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "missense")}
	if (GENETIC_CODE[[mut3]] == GENETIC_CODE[[codons[j]]]) {mutResult <- c(mutResult, "silent")}





j = j + 1

}

length(mutResult)
length(mutAlt)
length(refVec)
length(mutPos)


write.table(cbind(mutPos, refVec, mutAlt, mutResult), "lacZ_mutation_outcomes_result.txt", sep="\t", row.names=F)
