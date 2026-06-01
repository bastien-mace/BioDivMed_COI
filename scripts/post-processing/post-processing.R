# ---- 0. Load packages ----
library(tidyverse)
library(Biostrings) # handle DNA sequences
library(worrms) # extract WoRMS data
library(DECIPHER) # align DNA sequences
library(ape) # measure genetic distances
library(NeMO) # fit occupancy model
library(iNEXT) # estimate accumulation curves
library(rfishbase) # extract SeaLifeBase data
library(rentrez) # use NCBI 
library(taxize) # retrieve NCBI taxonomy
set_entrez_key('e0a14b0fc98a9389890653925881e0f04107') # NCBI API key
Sys.getenv('ENTREZ_KEY') 
# __________________________



# ---- 1. Import data ----
ASV_table <- read.csv('../../data/biodivmed_intra_metazoa_co1_dnoise_esv_table_rmtagjump.csv')
names(ASV_table) <- gsub('_MB0225D1', '', names(ASV_table))

MOTU_table <- read.csv('../../data/biodivmed_intra_metazoa_co1_swarm_mumu_motu_table.csv')
names(MOTU_table) <- gsub('_MB0225D1', '', names(MOTU_table))

## Write MOTUs and ASVs tables
write.csv(ASV_table, '../../data/ASV_table.csv', row.names = F)    
write.csv(MOTU_table, '../../data/MOTU_table.csv', row.names = F)
# ________________________



# ---- 2. Abundance filters and chimera removal ----
## 2-1. Remove ASVs following two criteria:
##  - ASVs representing less than 0.005% of total reads within each
##  - ASVs with less than 5 reads in total
ASV_filtering <- ASV_table %>%
  ### Only keep columns corresponding to samples
  pivot_longer(cols = starts_with('SPY'),
               names_to = 'sample',
               values_to = 'reads') %>%
  
  ### Calculate total number of reads by sample
  group_by(sample) %>%
  mutate(sample_total = sum(reads, na.rm = T),
         rel_abund = reads / sample_total) %>%
  
  ### First filter: Remove ESVs < 0.005% (0.00005) per sample
  mutate(reads = ifelse(rel_abund < 0.00005, 0, reads)) %>%
  ungroup() %>%
  
  ### Summarize by ASV after first filtering
  group_by(id) %>%
  summarise(total_reads_after_filter1 = sum(reads, na.rm = T), 
            .groups = 'drop') %>%
  
  ### Second filter: Remove ESVs with < 5 total reads
  filter(total_reads_after_filter1 >= 5)

### Recover ASV table with initial columns
ASV_table_filtered <- ASV_table %>% 
  filter(id %in% ASV_filtering$id)


## 2-2. Update MOTU table and remove chimeras
### Identify chimeras
MOTU_table_chimera <- read.csv('../../data/biodivmed_intra_metazoa_co1_table_motu.csv') %>% 
  filter(chimera != 'N')

### Filter MOTU table
MOTU_table_filtered <- MOTU_table %>% 
  filter(amplicon %in% ASV_table_filtered$MOTU,
         !amplicon %in% MOTU_table_chimera$amplicon) %>%
  dplyr::rename(seqid = amplicon) %>% 
  select(-c(CNEG01, best_identity, rank, taxid, scientific_name, species, genus, family, order, class))


## 2-3. Write MOTUs and ASVs tables
write.csv(ASV_table_filtered, '../../data/ASV_table_filtered.csv', row.names = F)    
write.csv(MOTU_table_filtered, '../../data/MOTU_table_filtered.csv', row.names = F)
# __________________________________________________



# ---- 3. Taxonomic assignment ----
## 3-1. Create a MOTU fasta file for taxonomic assignment with mkLTG
### Extract MOTUs' DNA sequences
dna_seqs <- DNAStringSet(MOTU_table_filtered$sequence)
names(dna_seqs) <- MOTU_table_filtered$seqid

### Write fasta file
writeXStringSet(dna_seqs, filepath = './taxonomy/MOTU_sequences.fasta', format = 'fasta')


## 3-2. Run taxonomic assignment with mkLTG on a Linux terminal
'
\bash
perl mkLTG/scripts/mkLTG.pl -in MOTU_sequences.fasta -blast_db COInr_filtered/COInr_filtered -taxonomy COInr_filtered/COInr_filtered_taxonomy.tsv -ltg_params mkLTG_params_modified_XT.tsv -outdir mkLTG_assignment/ -out_name mkLTG_assignment
'


## 3-3. Import assigned MOTUs
MOTU_table_assigned <- read_tsv('./taxonomy/mkLTG_assignment_ltg.tsv') %>% 
  select(-c(sequence)) %>% 
  inner_join(., MOTU_table_filtered) %>% 
  filter(!is.na(domain))
MOTU_table_assigned$kingdom <- gsub('kingdom_', '', MOTU_table_assigned$kingdom)
MOTU_table_assigned$kingdom <- gsub('phylum_', '', MOTU_table_assigned$kingdom)
MOTU_table_assigned$phylum <- gsub('phylum_', '', MOTU_table_assigned$phylum)
MOTU_table_assigned$class <- gsub('class_', '', MOTU_table_assigned$class)
MOTU_table_assigned$class <- gsub('order_', '', MOTU_table_assigned$class)
MOTU_table_assigned$order <- gsub('order_', '', MOTU_table_assigned$order)
MOTU_table_assigned$order <- gsub('family_', '', MOTU_table_assigned$order)
MOTU_table_assigned$family <- gsub('family_', '', MOTU_table_assigned$family)

### Restrict ASV table to assigned MOTUs
ASV_table_assigned <- ASV_table_filtered %>% 
  filter(MOTU %in% MOTU_table_assigned$seqid)


## 3-4. Write MOTUs and ASVs tables
write.csv(ASV_table_assigned, '../../data/ASV_table_assigned.csv', row.names = F)    
write.csv(MOTU_table_assigned, '../../data/MOTU_table_assigned.csv', row.names = F)
# _________________________________



# ---- 4. NUMTs correction ----
## 4-1. Prepare data
### 4-1-1. MOTU data
MOTU_table_to_correct <- MOTU_table_assigned %>% 
  dplyr::rename(id = seqid,
                COUNT = count,
                seq = sequence)
MOTU_table_to_correct$seq <- tolower(MOTU_table_to_correct$seq)

col_samples_MOTU <- c(15:32)

### 4-1-2. ASV data
ASV_table_to_correct <- ASV_table_assigned %>%
  dplyr::rename(COUNT = count,
                seq = sequence)
ASV_table_to_correct$seq <- tolower(ASV_table_to_correct$seq)

#### Some MOTU have an ID not corresponding to any ASV
no_representatives <- ASV_table_to_correct %>%
  filter(!MOTU %in% id)

#### Store the most abundant ASV within each of these MOTUs
correct_representative_ASVs <- no_representatives %>%
  group_by(MOTU) %>%
  slice_max(COUNT, n = 1, with_ties = F) %>% # keep the most abundant ASV
  ungroup() %>%
  select(correct_representative_ASV = id, MOTU)

### Update ASV table by correcting MOTU ID and sequence
ASV_table_to_correct <- ASV_table_to_correct %>%
  left_join(correct_representative_ASVs %>% # add a column with the corresponding MOTUs to its correct representative ASVs
              dplyr::rename(MOTU_rep = MOTU),
            by = c('id' = 'correct_representative_ASV')) %>%
  left_join(MOTU_table_to_correct %>% # add a column with the original MOTU sequence to the correct representative ASVs
      select(MOTU_rep = id, seq_MOTU = seq),
            by = 'MOTU_rep') %>%
  mutate(id  = if_else(!is.na(MOTU_rep), MOTU_rep, id), # if correct representative ASV: ASV ID become MOTU ID
         seq = if_else(!is.na(seq_MOTU), seq_MOTU, seq)) %>%  # if correct representative ASV: ASV sequence become MOTU sequence
  select(-MOTU_rep, -seq_MOTU)
ASV_table_to_correct$MOTU <- factor(ASV_table_to_correct$MOTU, levels = unique(ASV_table_to_correct$MOTU))

col_samples_ASV <- c(4:21)


## 4-2. Prepare NUMT tables
### 4-2-1. NUMT MOTU table
numts_MOTUs <- as.data.frame(matrix(NA, length(levels(ASV_table_to_correct$MOTU)), 17))

names(numts_MOTUs) <- c('MOTU_id', 'is_numt', 'is_metazoa_313',
                        'initial_ASV_richness', 'final_ASV_richness', 'ASVs_eliminated', 
                        'ASVs_no_stop', 'total_stops',  'ASVs_wrong_aa', 'ASVs_wrong_aa_no_stops')
numts_MOTUs$is_numt <- FALSE
numts_MOTUs$is_metazoa_313 <- FALSE
numts_MOTUs[, 4:8] <- 0

### 4-2-2. NUMT ASV table
numts_ASVs <- as.data.frame(matrix(NA, dim(ASV_table_to_correct)[1] ,10))

names(numts_ASVs) <- c('ASV_id', 'MOTU_id', 'is_numt', 'is_metazoa_313', 'n_stops',  'n_wrong_aa')
numts_ASVs$is_numt <- FALSE
numts_ASVs$is_metazoa_313 <- FALSE
numts_ASVs$n_stops <- 0


## 4-3. NUMT correction
for (i in 1:length(levels(ASV_table_to_correct$MOTU))){  # process MOTU by MOTU
  
  ### Extract ASVs belonging to the MOTU
  datas <- ASV_table_to_correct[which(ASV_table_to_correct$MOTU == levels(ASV_table_to_correct$MOTU)[i]), ] 
  
  ### Create a matrix to store ASVs's STOP codons
  stops <- matrix(NA, dim(datas)[1], 20) # 20 columns as up to 20 genetic codes are tested (see below)
  
  ### Transform ASVs' DNA sequences
  seq <- DNAStringSet(datas$seq)
  lenseq <- nchar(datas$seq)
  seq <- DNAStringSet(seq, start = 2, end = lenseq) # codon typically starts in nucleotide position 2 in the Leray's COI fragment
  
  ### Identify the representative ASV of the MOTU
  is_representative <- FALSE
  representative <- NULL
  for (j in 1:dim(datas)[1]){
    if (datas$seq[j] %in% MOTU_table_to_correct$seq){
      is_representative <- TRUE 
      representative <- j
    } 
  }
  if (is_representative == FALSE){
    message('ERROR: no representative ASV of MOTU ', MOTU_table_to_correct$id[i])
  }else{
    is_metazoa <- FALSE
  }
  
  ### Check if the MOTU is assigned to Metazoa kingdom
  if (!is.na(MOTU_table_to_correct$kingdom[which(MOTU_table_to_correct$seq == datas$seq[representative])])){
    if (MOTU_table_to_correct$kingdom[which(MOTU_table_to_correct$seq == datas$seq[representative])] == 'Metazoa'){
      is_metazoa <- TRUE 
    }
  }
  
  ### Translate ASVs' DNA sequences into amino-acid sequence
  if (is_metazoa){
    cod <- c(1, 2, 4, 5, 7, 11, 12, 15, 18) # for metazoans only genetic codes of metazoans are searched
  }else{
    cod <- c(1:20) # for non-metazoans the 20 genetic codes are searched
  }
  for (k in cod){
    code <- getGeneticCode(as.character(GENETIC_CODE_TABLE$id[k]))
    trans <- translate(seq, genetic.code = code) # get amino-acid sequence
    for (l in 1:dim(datas)[1]){
      aa <- strsplit(as.character(trans[l]), split = '')
      aa <- unlist(aa)
      nstops <- length(which(aa == '*')) # count number of STOP codons
      stops[l, k] <- nstops
    }
  }
  
  ### We retain codes with the minimum number of STOP codons
  best_code <- which(colSums(stops) == min(colSums(stops), na.rm = T))[1] # first code with minimum number of STOPs
  best_codes <- which(colSums(stops) == min(colSums(stops), na.rm = T)) # other codes with minimum number of STOPs
  
  ### For metazoans MOTU with 313 bp sequences, we check for wrong amino-acids in the five positions conserved according to Pentinsaari et al., 2016
  wrong_conserved_aa <- matrix(NA, dim(datas)[1], length(best_codes))
  if (is_metazoa & nchar(datas$seq[representative]) == 313){
    for (m in 1:length(best_codes)){
      code <- getGeneticCode(as.character(GENETIC_CODE_TABLE$id[best_codes[m]]))
      trans <- translate(seq, genetic.code = code) # get amino-acid sequence
      for (n in 1:dim(datas)[1]){
        aa <- strsplit(as.character(trans[n]), split = '')
        aa <- unlist(aa)
        wrong_aa <- 0
        if (aa[20] != 'H'){
          wrong_aa <- wrong_aa + 1
        } 
        if (aa[23] != 'G'){
          wrong_aa <- wrong_aa + 1
        } 
        if (aa[32] != 'N'){
          wrong_aa <- wrong_aa + 1
        } 
        if (aa[81] != 'D'){
          wrong_aa <- wrong_aa + 1
        } 
        if (aa[95] != 'G'){
          wrong_aa <- wrong_aa + 1
        } 
        wrong_conserved_aa[n, m] <- wrong_aa
      }
    }
  }
    
  ### We retain codes with the mininum number of wrong amino-acids over all MOTU ASVs'
  best_code_aa <- which(colSums(wrong_conserved_aa) == min(colSums(wrong_conserved_aa)))[1]  # first code with minimum number of wrong amino-acids
  
  ### Fill NUMT tables
  #### NUMT MOTU table
  numts_MOTUs$MOTU_id[i] <- datas$id[representative]
  numts_MOTUs$initial_ASV_richness[i] <- dim(datas)[1]
  numts_MOTUs$ASVs_no_stop[i] <- sum(stops[, best_code] == 0)
  numts_MOTUs$total_stops[i] <- sum(stops[, best_code])
  if (is_metazoa & nchar(datas$seq[representative]) == 313){
    numts_MOTUs$is_metazoa_313[i] <- TRUE
    numts_MOTUs$ASVs_wrong_aa[i] <- length(which(wrong_conserved_aa[, best_code_aa] > 0))
    numts_MOTUs$ASVs_wrong_aa_no_stops[i] <- length(which(stops[, best_code] == 0 & wrong_conserved_aa[, best_code_aa] > 0))
  } 
  #### ASVs with STOP codons and ASVs without STOP codons but with wrong amino-acids will be removed
  numts_MOTUs$final_ASV_richness[i] <- sum(numts_MOTUs$ASVs_no_stop[i], -numts_MOTUs$ASVs_wrong_aa_no_stops[i], na.rm = T)
  numts_MOTUs$ASVs_eliminated[i] <- numts_MOTUs$initial_ASV_richness[i] - numts_MOTUs$final_ASV_richness[i]
  if (numts_MOTUs$ASVs_eliminated[i] == numts_MOTUs$initial_ASV_richness[i]){
    numts_MOTUs$is_numt[i] <- TRUE # MOTU is a NUMT if all its ASVs are removed
  }
  
  #### NUMT ASV table
  positions <- which(ASV_table_to_correct$MOTU == levels(ASV_table_to_correct$MOTU)[i])
  flag <- vector('logical', length(positions))
  for (o in 1:length(positions)){
    pos <- positions[o]
    numts_ASVs$ASV_id[pos] <- ASV_table_to_correct$id[pos]
    numts_ASVs$MOTU_id[pos] <- datas$id[representative]
    numts_ASVs$n_stops[pos] <- stops[o, best_code]
    numts_ASVs$is_metazoa_313[pos] <- (is_metazoa & nchar(datas$seq[representative]) == 313)
    if (is_metazoa & nchar(datas$seq[representative]) == 313){
      numts_ASVs$n_wrong_aa[pos] <- wrong_conserved_aa[o, best_code_aa]
    }
    flag[o] <- sum(numts_ASVs$n_stops[pos], numts_ASVs$n_wrong_aa[pos], na.rm = T) > 0
    if (flag[o]){
      numts_ASVs$is_numt[pos] <- TRUE # ASV is a NUMT if sequence contains a STOP codon or a wrong amino-acid
    } 
  }
}


## 4-4. Remove NUMTs MOTUs
flag <- vector('logical', dim(numts_MOTUs)[1])
for (i in 1:dim(numts_MOTUs)[1]){
  MOTU <- numts_MOTUs[which(numts_MOTUs$MOTU_id == MOTU_table_to_correct$id[i]), ]
  if (MOTU$is_numt){
    flag[i] <- TRUE
  }
}
MOTU_table_corrected <- MOTU_table_to_correct
if (sum(flag) > 0){
  MOTU_table_corrected <- MOTU_table_corrected[-which(flag), ]
}


## 4-5. Remove NUMTs ASVs
flag <- vector('logical', dim(numts_ASVs)[1])
for (i in 1:dim(numts_ASVs)[1]){
  ASV <- numts_ASVs[which(numts_ASVs$ASV_id == ASV_table_to_correct$id[i]), ]
  if (ASV$is_numt){
    flag[i] <- TRUE
  }
}
ASV_table_corrected <- ASV_table_to_correct
if (sum(flag) > 0){
  ASV_table_corrected <- ASV_table_corrected[-which(flag), ]
}
ASV_table_corrected$MOTU <- factor(ASV_table_corrected$MOTU, levels = unique(ASV_table_corrected$MOTU))


## 4-6. Recalculate ASV abundance within MOTUs 
for (i in 1:length(levels(ASV_table_corrected$MOTU))){
  MOTU_table_corrected[which(MOTU_table_corrected$id == levels(ASV_table_corrected$MOTU)[i]), col_samples_MOTU] <- as.list(colSums(ASV_table_corrected[which(ASV_table_corrected$MOTU == levels(ASV_table_corrected$MOTU)[i]), col_samples_ASV]))
}
MOTU_table_corrected$COUNT<-rowSums(MOTU_table_corrected[, col_samples_MOTU])


## 4-7. Write MOTUs and ASVs tables
write.csv(ASV_table_corrected, '../../data/ASV_table_corrected.csv', row.names = F)    
write.csv(MOTU_table_corrected, '../../data/MOTU_table_corrected.csv', row.names = F)   
write.csv(numts_MOTUs, './numts/numts_MOTUs.csv', row.names = F)
write.csv(numts_ASVs, './numts/numts_ASVs.csv', row.names = F)
# _____________________________



# ---- 5. Marine assignation checking ----
## 5-1. Remove non marine kingdoms
kingdoms <- unique(na.omit(MOTU_table_marine$kingdom))

### Query WoRMS database
worms_query_kingdoms <- wm_records_names(kingdoms)

### Correct kingdoms
kingdoms_to_correct <- data.frame(initial_kingdom = character(),
                                  corrected_kingdom = character)
for (i in 1:length(kingdoms)){
  if (dim(worms_query_kingdoms[[i]])[1] == 1){
    if (worms_query_kingdoms[[i]]$rank != 'Kingdom'){
      correction <- data.frame(initial_kingdom = kingdoms[i],
                               corrected_kingdom = worms_query_kingdoms[[i]]$kingdom)
      kingdoms_to_correct <- rbind(kingdoms_to_correct,
                                   correction)
    }
    if (worms_query_kingdoms[[i]]$status == 'unaccepted'){
      correction <- data.frame(initial_kingdom = kingdoms[i],
                               corrected_kingdom = worms_query_kingdoms[[i]]$valid_name)
      kingdoms_to_correct <- rbind(kingdoms_to_correct,
                                   correction)
    }
  }else{
    print(kingdoms[i])
  }
}
kingdoms_to_correct <- kingdoms_to_correct[apply(kingdoms_to_correct, 1, function(x) length(unique(x)) > 1), ]
kingdoms_to_correct
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Phaeophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Bacillariophyta'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Rhodophyta'] <- 'Plantae'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Oomycota'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Viridiplantae'] <- 'Plantae'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Raphidophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Dinophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Cryptophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Cercozoa'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Haptophyta'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Eustigmatophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Discosea'] <- 'Protozoa'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Chrysophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Dictyochophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Pelagophyceae'] <- 'Chromista'
MOTU_table_marine$kingdom[MOTU_table_marine$kingdom == 'Synurophyceae'] <- 'Chromista'

### Store marine and non-marine kingdoms
kingdoms <- unique(na.omit(MOTU_table_marine$kingdom))
worms_query_kingdoms <- wm_records_names(kingdoms)
marine_kingdoms <- c()
non_marine_kingdoms <- c()
for (i in 1:length(kingdoms)){
  if (length(worms_query_kingdoms[[i]]) != 0){
    worms_query_kingdoms[[i]]$isMarine[is.na(worms_query_kingdoms[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_kingdoms[[i]]) == 0 || worms_query_kingdoms[[i]]$isMarine[1] != 1){
    non_marine_kingdoms <- c(non_marine_kingdoms, kingdoms[i])
  }else{
    marine_kingdoms <- c(marine_kingdoms, kingdoms[i])
  }
}
marine_kingdoms
non_marine_kingdoms
### NOTE: They are all marine !


## 5-2. Remove non marine phyla
phyla <- unique(na.omit(MOTU_table_marine$phylum))

### Query WoRMS database
worms_query_phyla <- wm_records_names(phyla)

### Correct phyla
phyla_to_correct <- data.frame(initial_phylum = character(),
                               corrected_phylum = character)
for (i in 1:length(phyla)){
  if (dim(worms_query_phyla[[i]])[1] == 1){
    if (worms_query_phyla[[i]]$rank != 'Phylum'){
      correction <- data.frame(initial_phylum = phyla[i],
                               corrected_phylum = worms_query_phyla[[i]]$phylum)
      phyla_to_correct <- rbind(phyla_to_correct,
                                correction)
    }
    if (worms_query_phyla[[i]]$status == 'unaccepted'){
      correction <- data.frame(initial_phylum = phyla[i],
                               corrected_phylum = worms_query_phyla[[i]]$valid_name)
      phyla_to_correct <- rbind(phyla_to_correct,
                                correction)
    }
  }else{
    print(phyla[i])
  }
}
phyla_to_correct <- phyla_to_correct[apply(phyla_to_correct, 1, function(x) length(unique(x)) > 1), ]
phyla_to_correct
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Phaeophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Raphidophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Dinophyceae'] <- 'Myzozoa'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Cryptophyceae'] <- 'Cryptophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Eustigmatophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Discosea'] <- 'Amoebozoa'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Chrysophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Dictyochophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Pelagophyceae'] <- 'Ochrophyta'
MOTU_table_marine$phylum[MOTU_table_marine$phylum == 'Synurophyceae'] <- 'Ochrophyta'

### Store marine and non-marine phyla
phyla <- unique(na.omit(MOTU_table_marine$phylum))
worms_query_phyla <- wm_records_names(phyla)
marine_phyla <- c()
non_marine_phyla <- c()
for (i in 1:length(phyla)){
  if (length(worms_query_phyla[[i]]) != 0){
    worms_query_phyla[[i]]$isMarine[is.na(worms_query_phyla[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_phyla[[i]]) == 0 || worms_query_phyla[[i]]$isMarine[1] != 1){
    non_marine_phyla <- c(non_marine_phyla, phyla[i])
  }else{
    marine_phyla <- c(marine_phyla, phyla[i])
  }
}
marine_phyla
non_marine_phyla
### NOTE: They are all marine !


## 5-3. Remove non marine classes
classes <- unique(na.omit(MOTU_table_marine$class))

### Query WoRMS database
worms_query_classes <- wm_records_names(classes)

### Correct classes
classes_to_correct <- data.frame(initial_class = character(),
                                 corrected_class = character)
for (i in 1:length(classes)){
  if (dim(worms_query_classes[[i]])[1] == 1){
    if (worms_query_classes[[i]]$rank != 'Class'){
      correction <- data.frame(initial_class = classes[i],
                               corrected_class = worms_query_classes[[i]]$class)
      classes_to_correct <- rbind(classes_to_correct,
                                  correction)
    }
    if (worms_query_classes[[i]]$status == 'unaccepted'){
      correction <- data.frame(initial_class = classes[i],
                               corrected_class = worms_query_classes[[i]]$valid_name)
      classes_to_correct <- rbind(classes_to_correct,
                                  correction)
    }
  }else{
    print(classes[i])
  }
}
classes_to_correct <- classes_to_correct[apply(classes_to_correct, 1, function(x) length(unique(x)) > 1), ]
classes_to_correct
MOTU_table_marine$class[MOTU_table_marine$class == 'Anthozoa' & MOTU_table_marine$order == 'Actiniaria'] <- 'Hexacorallia'
MOTU_table_marine$class[MOTU_table_marine$class == 'Anthozoa' & MOTU_table_marine$order == 'Scleractinia'] <- 'Hexacorallia'
MOTU_table_marine$class[MOTU_table_marine$class == 'Anthozoa' & MOTU_table_marine$order == 'Zoantharia'] <- 'Hexacorallia'
MOTU_table_marine$class[MOTU_table_marine$class == 'Anthozoa' & MOTU_table_marine$order == 'Malacalcyonacea'] <- 'Octocorallia'
MOTU_table_marine$class[MOTU_table_marine$class == 'Enopla'] <- 'Hoplonemertea'
MOTU_table_marine$class[MOTU_table_marine$class == 'Coscinodiscophyceae'] <- 'Bacillariophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Mediophyceae'] <- 'Bacillariophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Hexanauplia'] <- 'Copepoda'
MOTU_table_marine$class[MOTU_table_marine$class == 'Pythiales'] <- 'Peronosporea'
MOTU_table_marine$class[MOTU_table_marine$class == 'Saprolegniales'] <- 'Peronosporea'
MOTU_table_marine$class[MOTU_table_marine$class == 'Isochrysidales'] <- 'Coccolithophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Prymnesiales'] <- 'Coccolithophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Peronosporales'] <- 'Peronosporea'
MOTU_table_marine$class[MOTU_table_marine$class == 'Flabellinia'] <- 'Discosea'
MOTU_table_marine$class[MOTU_table_marine$class == 'Fragilariophyceae'] <- 'Bacillariophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Synurophyceae'] <- 'Chrysophyceae'
MOTU_table_marine$class[MOTU_table_marine$class == 'Pseudoscourfieldiophyceae'] <- 'Pyramimonadophyceae'

### Store marine and non-marine classes
classes <- unique(na.omit(MOTU_table_marine$class))
worms_query_classes <- wm_records_names(classes)
marine_classes <- c()
non_marine_classes <- c()
for (i in 1:length(classes)){
  if (length(worms_query_classes[[i]]) != 0){
    worms_query_classes[[i]]$isMarine[is.na(worms_query_classes[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_classes[[i]]) == 0 || worms_query_classes[[i]]$isMarine[1] != 1){
    non_marine_classes <- c(non_marine_classes, classes[i])
  }else{
    marine_classes <- c(marine_classes, classes[i])
  }
}
marine_classes
non_marine_classes
### NOTE: They are all marine !


## 5-4. Remove non marine orders
orders <- unique(na.omit(MOTU_table_marine$order))

### Query WoRMS database
worms_query_orders <- wm_records_names(orders)

### Correct orders
orders_to_correct <- data.frame(initial_order = character(),
                                corrected_order = character)
for (i in 1:length(orders)){
  if (dim(worms_query_orders[[i]])[1] == 1){
    if (worms_query_orders[[i]]$rank != 'order'){
      correction <- data.frame(initial_order = orders[i],
                               corrected_order = worms_query_orders[[i]]$order)
      orders_to_correct <- rbind(orders_to_correct,
                                 correction)
    }
    if (worms_query_orders[[i]]$status == 'unaccepted'){
      correction <- data.frame(initial_order = orders[i],
                               corrected_order = worms_query_orders[[i]]$valid_name)
      orders_to_correct <- rbind(orders_to_correct,
                                 correction)
    }
  }else{
    print(orders[i])
  }
}
orders_to_correct <- orders_to_correct[apply(orders_to_correct, 1, function(x) length(unique(x)) > 1), ]
orders_to_correct
MOTU_table_marine$order[MOTU_table_marine$order == 'Diplostraca' & MOTU_table_marine$family == 'Podonidae'] <- 'Onychopoda'
MOTU_table_marine$order[MOTU_table_marine$order == 'Diplostraca' & MOTU_table_marine$family == 'Sididae'] <- 'Ctenopoda'
MOTU_table_marine$order[MOTU_table_marine$order == 'Cerithiidae'] <- 'Caenogastropoda'
MOTU_table_marine$order[MOTU_table_marine$order == 'Poecilostomatoida'] <- 'Cyclopoida'
MOTU_table_marine$order[MOTU_table_marine$order == 'Bigelowiella'] <- 'Chlorarachniales'
MOTU_table_marine$order[MOTU_table_marine$order == 'Plakobranchidae'] <- 'Sacoglossa'
MOTU_table_marine$order[MOTU_table_marine$order == 'Limapontiidae'] <- 'Sacoglossa'
MOTU_table_marine$order[MOTU_table_marine$order == 'Chaetocerotales'] <- 'Chaetocerotanae'
MOTU_table_marine$order[MOTU_table_marine$order == 'Valvatidae'] <- 'Heterobranchia'
MOTU_table_marine$order[MOTU_table_marine$order == 'Thraciidae'] <- 'Anomalodesmata'
MOTU_table_marine$order[MOTU_table_marine$order == 'Orbiniidae'] <- 'Scolecida'
MOTU_table_marine$order[MOTU_table_marine$order == 'Paraonidae'] <- 'Scolecida'
MOTU_table_marine$order[MOTU_table_marine$order == 'Moronidae'] <- 'Perciformes'
MOTU_table_marine$order[MOTU_table_marine$order == 'Maldanidae'] <- 'Scolecida'
MOTU_table_marine$order[MOTU_table_marine$order == 'Chaetopteridae'] <- 'Sedentaria'
MOTU_table_marine$order[MOTU_table_marine$order == 'Stephanodiscales'] <- 'Thalassiosiranae'
MOTU_table_marine$order[MOTU_table_marine$order == 'Cyathomonadacea'] <- 'Cryptophyceae'

### Store marine and non-marine orders
orders <- unique(na.omit(MOTU_table_marine$order))
worms_query_orders <- wm_records_names(orders)
marine_orders <- c()
non_marine_orders <- c()
for (i in 1:length(orders)){
  if (length(worms_query_orders[[i]]) != 0){
    worms_query_orders[[i]]$isMarine[is.na(worms_query_orders[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_orders[[i]]) == 0 || worms_query_orders[[i]]$isMarine[1] != 1){
    non_marine_orders <- c(non_marine_orders, orders[i])
  }else{
    marine_orders <- c(marine_orders, orders[i])
  }
}
marine_orders
non_marine_orders
### NOTE: One 'non-marine' orders is marine !
non_marine_orders <-  non_marine_orders[!non_marine_orders %in% c('Spariformes')]

### Remove non-marine orders
MOTU_table_marine <- MOTU_table_marine %>% 
  filter(!order %in% non_marine_orders)
ASV_table_marine <- ASV_table_marine %>% 
  filter(MOTU %in% MOTU_table_marine$id)


## 5-5. Remove non marine families
families <- unique(na.omit(MOTU_table_marine$family))

### Query WoRMS database
worms_query_families <- wm_records_names(families)

### Correct families
families_to_correct <- data.frame(initial_family = character(),
                                  corrected_family = character)
for (i in 1:length(families)){
  if (dim(worms_query_families[[i]])[1] == 1){
    if (worms_query_families[[i]]$rank != 'family'){
      correction <- data.frame(initial_family = families[i],
                               corrected_family = worms_query_families[[i]]$family)
      families_to_correct <- rbind(families_to_correct,
                                   correction)
    }
    if (worms_query_families[[i]]$status == 'unaccepted'){
      correction <- data.frame(initial_family = families[i],
                               corrected_family = worms_query_families[[i]]$valid_name)
      families_to_correct <- rbind(families_to_correct,
                                   correction)
    }
  }else{
    print(families[i])
  }
}
families_to_correct <- families_to_correct[apply(families_to_correct, 1, function(x) length(unique(x)) > 1), ]
families_to_correct
MOTU_table_marine$family[MOTU_table_marine$family == 'Bigelowiella'] <- 'Chlorarachniaceae'
MOTU_table_marine$family[MOTU_table_marine$family == 'Euterpinidae'] <- 'Tachidiidae'
MOTU_table_marine$family[MOTU_table_marine$family == 'Quasitetrastemma'] <- 'Tetrastemmatidae'
MOTU_table_marine$family[MOTU_table_marine$family == 'Pseudochattonella'] <- 'Sacoglossa'
MOTU_table_marine$family[MOTU_table_marine$family == 'Skeletonemataceae'] <- 'Skeletonemaceae'
MOTU_table_marine$family[MOTU_table_marine$family == 'Clytiidae'] <- 'Campanulariidae'

### Store marine and non-marine families
families <- unique(na.omit(MOTU_table_marine$family))
worms_query_families <- wm_records_names(families)
marine_families <- c()
non_marine_families <- c()
for (i in 1:length(families)){
  if (length(worms_query_families[[i]]) != 0){
    worms_query_families[[i]]$isMarine[is.na(worms_query_families[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_families[[i]]) == 0 || worms_query_families[[i]]$isMarine[1] != 1){
    non_marine_families <- c(non_marine_families, families[i])
  }else{
    marine_families <- c(marine_families, families[i])
  }
}
marine_families
non_marine_families

### Remove non-marine families
MOTU_table_marine <- MOTU_table_marine %>% 
  filter(!family %in% non_marine_families)
ASV_table_marine <- ASV_table_marine %>% 
  filter(MOTU %in% MOTU_table_marine$id)


## 5-6. Remove non marine genera
genera <- unique(na.omit(MOTU_table_marine$genus))

### Query WoRMS database
worms_query_genera <- wm_records_names(genera)

### Store marine and non-marine genera
marine_genera <- c()
non_marine_genera <- c()
for (i in 1:length(genera)){
  if (length(worms_query_genera[[i]]) != 0){
    worms_query_genera[[i]]$isMarine[is.na(worms_query_genera[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_genera[[i]]) == 0 || worms_query_genera[[i]]$isMarine[1] != 1){
    non_marine_genera <- c(non_marine_genera, genera[i])
  }else{
    marine_genera <- c(marine_genera, genera[i])
  }
}
marine_genera
non_marine_genera
### NOTE: One 'non-marine' genus is marine !
non_marine_genera <-  non_marine_genera[!non_marine_genera %in% c('Stephanodiscus')]

### Remove non-marine genera
MOTU_table_marine <- MOTU_table_marine %>% 
  filter(!genus %in% non_marine_genera)
ASV_table_marine <- ASV_table_marine %>% 
  filter(MOTU %in% MOTU_table_marine$id)


## 5-7. Remove non marine species
species <- unique(na.omit(MOTU_table_marine$species))

### Query WoRMS database
worms_query_species <- wm_records_names(species)

### Store marine and non-marine species
marine_species <- c()
non_marine_species <- c()
for (i in 1:length(species)){
  if (length(worms_query_species[[i]]) != 0){
    worms_query_species[[i]]$isMarine[is.na(worms_query_species[[i]]$isMarine)] <- 0
  }
  if (length(worms_query_species[[i]]) == 0 || worms_query_species[[i]]$isMarine[1] != 1){
    non_marine_species <- c(non_marine_species, species[i])
  }else{
    marine_species <- c(marine_species, species[i])
  }
}
marine_species
non_marine_species
### NOTE: One 'non-marine' species is marine !
non_marine_species <-  non_marine_species[!non_marine_species %in% c('Lithodesmium intricatum')]

### Remove non-marine species
MOTU_table_marine <- MOTU_table_marine %>% 
  filter(!species %in% non_marine_species)
ASV_table_marine <- ASV_table_marine %>% 
  filter(MOTU %in% MOTU_table_marine$id)


## 5-8. Write MOTUs and ASVs tables
write.csv(ASV_table_marine, '../../data/ASV_table_marine.csv', row.names = F)    
write.csv(MOTU_table_marine, '../../data/MOTU_table_marine.csv', row.names = F)
# ________________________________________



# ---- 6. Barcoding gap analysis ----
## 6-1. Inter-MOTU level
### Extract MOTU and align sequences
MOTU_sequences <- DNAStringSet(MOTU_table_marine$seq)
names(MOTU_sequences) <- MOTU_table_marine$id
MOTU_sequences_aligned <- AlignSeqs(MOTU_sequences)

### Calculate nucleotide distance
inter_MOTU_dist <- dist.dna(as.DNAbin(MOTU_sequences_aligned),
                            model = 'raw')

### Summary
cat('Mean inter-MOTU distance =', mean(inter_MOTU_dist), '\n',
    '95% CI = [', quantile(inter_MOTU_dist, 0.025), '-', quantile(inter_MOTU_dist, 0.975), ']\n',
    'Number of MOTU pairs =', length(inter_MOTU_dist), '\n')


## 6-2. Intra-MOTU level
### Identify MOTU with variability
MOTU_with_different_ASVs <- ASV_table_marine %>%
  count(MOTU, name = 'n_ASV') %>% 
  dplyr::rename(id = MOTU) %>% 
  filter(n_ASV > 1)

### Extract corresponding ASVs
ASV_restricted <- ASV_table_marine %>% 
  filter(MOTU %in% MOTU_with_different_ASVs$id)

### Measure intra-MOTU variability
nb <- nrow(MOTU_with_different_ASVs)
intra_MOTU_dist <- vector('list', nb)
names(intra_MOTU_dist) <- MOTU_with_different_ASVs$id
for (i in 1:nb){
  MOTU_i <- MOTU_with_different_ASVs$id[i]
  message('MOTU ', i, '/', nb, ' : ', MOTU_i)
  
  ### Extract MOTU ASVs and align ASV sequences
  intra_MOTU_i <- ASV_table_marine %>%
    filter(MOTU == MOTU_i)
  ASV_sequences <- DNAStringSet(intra_MOTU_i$seq)
  names(ASV_sequences) <- ASV_sequences$id
  
  ### Calculate nucleotide distance
  ASV_sequences_aligned <- AlignSeqs(ASV_sequences)
  intra_MOTU_dist[[i]] <- dist.dna(as.DNAbin(ASV_sequences_aligned),
                                   model = 'raw')
}

### Extract all values
all_intra_MOTU_dist <- unlist(intra_MOTU_dist)

### Summary
cat('Mean intra-MOTU distance =', mean(all_intra_MOTU_dist), '\n',
    '95% CI = [', quantile(all_intra_MOTU_dist, 0.025), '-', quantile(all_intra_MOTU_dist, 0.975), ']\n',
    'Number of ASV pairs =', length(all_intra_MOTU_dist), '\n')


## 6-3. Check for barcoding gap
distance_df <- data.frame(distance = c(as.vector(inter_MOTU_dist),
                                       as.vector(all_intra_MOTU_dist)),
                          type = c(rep('Inter-MOTU', length(inter_MOTU_dist)),
                                   rep('Intra-MOTU', length(all_intra_MOTU_dist))))

### Plot histograms
barcoding_gap_before <- ggplot(distance_df, aes(x = distance, fill = type)) +
  geom_histogram(alpha = 0.5, position = 'identity', bins = 40) +
  scale_x_continuous(breaks = seq(0, max(distance_df$distance), by = 0.1)) +
  labs(x = 'Pairwise genetic distance (p-distance)',
       y = 'Frequency',
       fill = '') +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.3, 0.9),
        legend.justification = c(0.5, 1))
barcoding_gap_before

### Save in pdf
ggsave(barcoding_gap_before,
       file = './barcoding_gap/barcoding_gap_before.pdf',
       height = 10, width = 15, unit = 'cm')
### NOTE: Some intra-MOTU distance overlap inter-MOTU distances !
### Those hypervariable MOTUs should probably be splitted into different MOTUs


## 6-4. Deal with hypervariable MOTUs
### We identify MOTUs with genetic distance above 0.15
hypervariation <- c()
for(i in 1:length(intra_MOTU_dist)){
  if(sum(unlist(intra_MOTU_dist[i]) > 0.15) > 1){
    hypervariation <- c(hypervariation, i)
  }
}
hypervariable_MOTU <- MOTU_table_marine %>% 
  filter(id %in% names(intra_MOTU_dist)[hypervariation])
### NOTE: 2 MOTUs are concerned

### Plot their distance profiles
pdf('./barcoding_gap/hypervariable_MOTUs.pdf',
     height = 10, width = 10)
par(mfrow = c(2,1))
for(i in hypervariable_MOTU$id){
  x <- as.numeric(intra_MOTU_dist[[i]])
  hist(x, main = i,  breaks = 20,
       xlab = 'Distance', xlim = c(0, max(all_intra_MOTU_dist)))
}
### MOTUs 61eae644abf0141da7db5cd5597786614c7ef9ce and fec200bfb97c8891c446cf82d04c3943669a5d2b
### show bimodal distributions with one peak around 0.15-0.2 !
### They are likely composed of different MOTUs, not differentiated enough for the clustering algorithm used !
### We will split them as this may better reflect biological reality.
dev.off()


## 6-5. Splitting hypervariable MOTUs
### 6-5-1. MOTU 61eae644abf0141da7db5cd5597786614c7ef9ce
intra_MOTU_dist_61e <- intra_MOTU_dist$`61eae644abf0141da7db5cd5597786614c7ef9ce`

#### 6-5-1-1. Identify adequate clusters
##### Represent a dendogram based on sequence similarity
pdf('./barcoding_gap/61e_MOTU_n_cluster_2.pdf',
     height = 10, width = 10)
layout(matrix(c(1, 1,
                2, 3),
              nrow = 2, byrow = TRUE))
hc_61e <- hclust(intra_MOTU_dist_61e, method = 'average')
plot(hc_61e, labels = F,
     ylab = 'Distance', xlab = '', sub = '',
     main = 'Dendogram of ASV sequences\nfor MOTU 61eae644abf0141da7db5cd5597786614c7ef9ce')

##### We try 2 clusters
rect.hclust(hc_61e, k = 2)
clusters_61e_n2 <- cutree(hc_61e, k = 2)
table(clusters_61e_n2)

##### Calculate intra-cluster genetic distances
new_intra_MOTU_dist_61e <- as.matrix(intra_MOTU_dist_61e)
idx1_61e_n2 <- which(clusters_61e_n2 == 1)
idx2_61e_n2 <- which(clusters_61e_n2 == 2)

###### Cluster 1
intra_matrix_cluster_1_61e_n2 <- new_intra_MOTU_dist_61e[idx1_61e_n2, idx1_61e_n2]
intra_dist_cluster_1_61e_n2 <- intra_matrix_cluster_1_61e_n2[upper.tri(intra_matrix_cluster_1_61e_n2)]
summary(intra_dist_cluster_1_61e_n2)
hist(intra_dist_cluster_1_61e_n2, main = 'Cluster 1', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))
###### NOTE: There is still high values

###### Cluster 2
intra_matrix_cluster_2_61e_n2 <- new_intra_MOTU_dist_61e[idx2_61e_n2, idx2_61e_n2]
intra_dist_cluster_2_61e_n2 <- intra_matrix_cluster_2_61e_n2[upper.tri(intra_matrix_cluster_2_61e_n2)]
summary(intra_dist_cluster_2_61e_n2)
hist(intra_dist_cluster_2_61e_n2, main = 'Cluster 2', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))
dev.off()

##### We try 3 clusters
pdf('./barcoding_gap/61e_MOTU_n_cluster_3.pdf',
     height = 10, width = 10)
layout(matrix(c(1, 1, 1,
                2, 3, 4),
              nrow = 2, byrow = TRUE))
plot(hc_61e, labels = F,
     ylab = 'Distance', xlab = '', sub = '',
     main = 'Dendogram of ASV sequences\nfor MOTU 61eae644abf0141da7db5cd5597786614c7ef9ce')
rect.hclust(hc_61e, k = 3)
clusters_61e_n3 <- cutree(hc_61e, k = 3)
table(clusters_61e_n3)

##### Calculate intra-cluster genetic distances
new_intra_MOTU_dist_61e <- as.matrix(intra_MOTU_dist_61e)
idx1_61e_n3 <- which(clusters_61e_n3 == 1)
idx2_61e_n3 <- which(clusters_61e_n3 == 2)
idx3_61e_n3 <- which(clusters_61e_n3 == 3)

###### Cluster 1
intra_matrix_cluster_1_61e_n3 <- new_intra_MOTU_dist_61e[idx1_61e_n3, idx1_61e_n3]
intra_dist_cluster_1_61e_n3 <- intra_matrix_cluster_1_61e_n3[upper.tri(intra_matrix_cluster_1_61e_n3)]
summary(intra_dist_cluster_1_61e_n3)
hist(intra_dist_cluster_1_61e_n3, main = 'Cluster 1', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))

###### Cluster 2
intra_matrix_cluster_2_61e_n3 <- new_intra_MOTU_dist_61e[idx2_61e_n3, idx2_61e_n3]
intra_dist_cluster_2_61e_n3 <- intra_matrix_cluster_2_61e_n3[upper.tri(intra_matrix_cluster_2_61e_n3)]
summary(intra_dist_cluster_2_61e_n3)
hist(intra_dist_cluster_2_61e_n3, main = 'Cluster 2', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))

###### Cluster 3
intra_matrix_cluster_3_61e_n3 <- new_intra_MOTU_dist_61e[idx3_61e_n3, idx3_61e_n3]
intra_dist_cluster_3_61e_n3 <- intra_matrix_cluster_3_61e_n3[upper.tri(intra_matrix_cluster_3_61e_n3)]
summary(intra_dist_cluster_3_61e_n3)
hist(intra_dist_cluster_3_61e_n3, main = 'Cluster 3', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))
dev.off()

#### 6-5-1-2. Update ASV and MOTU tables
MOTU_61e <- '61eae644abf0141da7db5cd5597786614c7ef9ce'

##### Associate each ASV to its cluster
ASV_61e <- ASV_table_marine %>%
  filter(MOTU == MOTU_61e) %>% 
  mutate(cluster = clusters_61e_n3)

##### Extract IDs and sequences of new MOTUs
sample_cols_61e <- grep('^SPY', names(ASV_61e), value = T)

new_MOTUs_61e <- ASV_61e %>%
  group_by(cluster) %>%
  summarise(new_MOTU_id  = id[which.max(COUNT)],
    new_MOTU_seq = seq[which.max(COUNT)],
    COUNT = sum(COUNT),
    across(all_of(sample_cols_61e), sum), .groups = 'drop') %>%
  mutate(old_MOTU = MOTU_61e)

##### Identify to which new MOTU ASVs belong
map_ASV_to_new_MOTU_61e <- ASV_61e %>%
  select(id, cluster) %>%
  left_join(new_MOTUs_61e %>% select(cluster, new_MOTU_id),
            by = 'cluster')

##### Update ASV table
ASV_table_updated <- ASV_table_marine %>%
  left_join(map_ASV_to_new_MOTU_61e %>% select(id, new_MOTU_id),
            by = 'id') %>%
  mutate(MOTU = ifelse(MOTU == MOTU_61e,
                       new_MOTU_id,
                       MOTU)) %>%
  select(-new_MOTU_id)

##### Update MOTU table
MOTU_table_updated <- MOTU_table_marine %>%
  filter(id != MOTU_61e)
MOTU_table_updated <- bind_rows(MOTU_table_updated,
                                       new_MOTUs_61e %>%
                                         mutate(id = new_MOTU_id,
                                                seq = new_MOTU_seq)) %>% 
  select(-c(cluster, new_MOTU_id, new_MOTU_seq, old_MOTU))

### 6-5-2. MOTU fec200bfb97c8891c446cf82d04c3943669a5d2b
intra_MOTU_dist_fec <- intra_MOTU_dist$`fec200bfb97c8891c446cf82d04c3943669a5d2b`

#### 6-5-2-1. Identify adequate clusters
##### Represent a dendogram based on sequence similarity
pdf('./barcoding_gap/fec_MOTU_n_cluster_2.pdf',
    height = 10, width = 10)
layout(matrix(c(1, 1,
                2, 3),
              nrow = 2, byrow = TRUE))
hc_fec <- hclust(intra_MOTU_dist_fec, method = 'average')
plot(hc_fec, labels = F,
     ylab = 'Distance', xlab = '', sub = '',
     main = 'Dendogram of ASV sequences\nfor MOTU fec200bfb97c8891c446cf82d04c3943669a5d2b')

##### We try 2 clusters
rect.hclust(hc_fec, k = 2)
clusters_fec_n2 <- cutree(hc_fec, k = 2)
table(clusters_fec_n2)

##### Calculate intra-cluster genetic distances
new_intra_MOTU_dist_fec <- as.matrix(intra_MOTU_dist_fec)
idx1_fec_n2 <- which(clusters_fec_n2 == 1)
idx2_fec_n2 <- which(clusters_fec_n2 == 2)

###### Cluster 1
intra_matrix_cluster_1_fec_n2 <- new_intra_MOTU_dist_fec[idx1_fec_n2, idx1_fec_n2]
intra_dist_cluster_1_fec_n2 <- intra_matrix_cluster_1_fec_n2[upper.tri(intra_matrix_cluster_1_fec_n2)]
summary(intra_dist_cluster_1_fec_n2)
hist(intra_dist_cluster_1_fec_n2, main = 'Cluster 1', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))
###### NOTE: There is still high values

###### Cluster 2
intra_matrix_cluster_2_fec_n2 <- new_intra_MOTU_dist_fec[idx2_fec_n2, idx2_fec_n2]
intra_dist_cluster_2_fec_n2 <- intra_matrix_cluster_2_fec_n2[upper.tri(intra_matrix_cluster_2_fec_n2)]
summary(intra_dist_cluster_2_fec_n2)
hist(intra_dist_cluster_2_fec_n2, main = 'Cluster 2', xlab = 'Intra-MOTU distance', xlim = c(0, max(all_intra_MOTU_dist)))
dev.off()

#### 6-5-2-2. Update ASV and MOTU tables
MOTU_fec <- 'fec200bfb97c8891c446cf82d04c3943669a5d2b'

##### Associate each ASV to its cluster
ASV_fec <- ASV_table_marine %>%
  filter(MOTU == MOTU_fec) %>% 
  mutate(cluster = clusters_fec_n2)

##### Extract IDs and sequences of new MOTUs
sample_cols_fec <- grep('^SPY', names(ASV_fec), value = T)

new_MOTUs_fec <- ASV_fec %>%
  group_by(cluster) %>%
  summarise(new_MOTU_id  = id[which.max(COUNT)],
            new_MOTU_seq = seq[which.max(COUNT)],
            COUNT = sum(COUNT),
            across(all_of(sample_cols_fec), sum), .groups = 'drop') %>%
  mutate(old_MOTU = MOTU_fec)

##### Identify to which new MOTU ASVs belong
map_ASV_to_new_MOTU_fec <- ASV_fec %>%
  select(id, cluster) %>%
  left_join(new_MOTUs_fec %>% select(cluster, new_MOTU_id),
            by = 'cluster')

##### Update ASV table
ASV_table_updated <- ASV_table_updated %>%
  left_join(map_ASV_to_new_MOTU_fec %>% select(id, new_MOTU_id),
            by = 'id') %>%
  mutate(MOTU = ifelse(MOTU == MOTU_fec,
                       new_MOTU_id,
                       MOTU)) %>%
  select(-new_MOTU_id)

##### Update MOTU table
MOTU_table_updated <- MOTU_table_updated %>%
  filter(id != MOTU_fec)
MOTU_table_updated <- bind_rows(MOTU_table_updated,
                                       new_MOTUs_fec %>%
                                         mutate(id = new_MOTU_id,
                                                seq = new_MOTU_seq)) %>% 
  select(-c(cluster, new_MOTU_id, new_MOTU_seq, old_MOTU))


## 6-6. Rerun taxonomic assignment for new MOTUs
### 6-6-1. Create a MOTU fasta file for taxonomic assignment with mkLTG
#### Extract MOTUs' DNA sequences
new_MOTU_seqs <- DNAStringSet(MOTU_table_updated$seq[is.na(MOTU_table_updated$ltg_taxid)])
names(new_MOTU_seqs) <- MOTU_table_updated$id[is.na(MOTU_table_updated$ltg_taxid)]

#### Write fasta file
# writeXStringSet(new_MOTU_seqs, filepath = './taxonomy/new_MOTU_sequences.fasta', format = 'fasta')

### 6-6-2. Run taxonomic assignment with mkLTG
'
\bash
perl mkLTG/scripts/mkLTG.pl -in new_MOTU_sequences.fasta -blast_db COInr_filtered/COInr_filtered -taxonomy COInr_filtered/COInr_filtered_taxonomy.tsv -ltg_params mkLTG_params_modified_XT.tsv -outdir new_mkLTG_assignment/ -out_name new_mkLTG_assignment
'

### 6-6-3. Import assigned MOTUs
new_assignment <- read_tsv('./taxonomy/new_mkLTG_assignment_ltg.tsv') %>% 
  dplyr::rename(id = seqid,
                seq = sequence)
new_assignment$kingdom <- gsub('kingdom_', '', new_assignment$kingdom)

#### Correct taxonomy
new_assignment$kingdom[new_assignment$kingdom == 'Bacillariophyta'] <- 'Chromista'
new_assignment$class[new_assignment$class == 'Coscinodiscophyceae'] <- 'Bacillariophyceae'
new_assignment$family[new_assignment$family == 'Skeletonemataceae'] <- 'Skeletonemaceae'

#### Add the new MOTUs
MOTU_table_updated <- left_join(MOTU_table_updated, new_assignment %>% 
                                  select(id, pid, ltg_taxid, ltg_name, ltg_rank, ltg_rank_index, domain, kingdom, phylum, class, order, family, genus, species),
                                by = 'id',
                                suffix = c('', '.new')) %>% 
  mutate(pid = coalesce(pid.new, pid),
         ltg_taxid = coalesce(ltg_taxid.new, ltg_taxid),
         ltg_name = coalesce(ltg_name.new, ltg_name),
         ltg_rank = coalesce(ltg_rank.new, ltg_rank),
         ltg_rank_index = coalesce(ltg_rank_index.new, ltg_rank_index),
         domain = coalesce(domain.new, domain),
         kingdom = coalesce(kingdom.new, kingdom),
         phylum = coalesce(phylum.new, phylum),
         class = coalesce(class.new, class),
         order = coalesce(order.new, order),
         family = coalesce(family.new, family),
         genus = coalesce(genus.new, genus),
         species = coalesce(species.new, species)) %>% 
  select(-ends_with('.new'))


## 6-7. Recalculate inter-MOTU and intra-MOTU distances 
### 6-7-1. Inter-MOTU level
#### Extract MOTU and align sequences
MOTU_sequences_new <- DNAStringSet(MOTU_table_updated$seq)
names(MOTU_sequences_new) <- MOTU_table_updated$id
MOTU_sequences_aligned_new <- AlignSeqs(MOTU_sequences_new)

#### Calculate nucleotide distance
inter_MOTU_dist_new <- dist.dna(as.DNAbin(MOTU_sequences_aligned_new),
                                model = 'raw')

#### Summary
cat('Mean inter-MOTU distance =', mean(inter_MOTU_dist_new), '\n',
    '95% CI = [', quantile(inter_MOTU_dist_new, 0.025), '-', quantile(inter_MOTU_dist_new, 0.975), ']\n',
    'Number of MOTU pairs =', length(inter_MOTU_dist_new), '\n')

### 6-7-2. Intra-MOTU level
#### Identify MOTU with variability
MOTU_with_different_ASVs_new <- ASV_table_updated %>%
  count(MOTU, name = 'n_ASV') %>% 
  dplyr::rename(id = MOTU) %>% 
  filter(n_ASV > 1)

#### Extract corresponding ASVs
ASV_restricted_new <- ASV_table_updated %>% 
  filter(MOTU %in% MOTU_with_different_ASVs_new$id)

#### Measure intra-MOTU variability
nb_new <- nrow(MOTU_with_different_ASVs_new)
intra_MOTU_dist_new <- vector('list', nb_new)
names(intra_MOTU_dist_new) <- MOTU_with_different_ASVs_new$id
for (i in 1:nb_new){
  MOTU_i <- MOTU_with_different_ASVs_new$id[i]
  message('MOTU ', i, '/', nb_new, ' : ', MOTU_i)
  
  #### Extract MOTU ASVs and align ASV sequences
  intra_MOTU_i <- ASV_table_updated %>%
    filter(MOTU == MOTU_i)
  ASV_sequences <- DNAStringSet(intra_MOTU_i$seq)
  names(ASV_sequences) <- ASV_sequences$id
  
  #### Calculate nucleotide distance
  ASV_sequences_aligned <- AlignSeqs(ASV_sequences)
  intra_MOTU_dist_new[[i]] <- dist.dna(as.DNAbin(ASV_sequences_aligned),
                                       model = 'raw')
}

#### Extract all values
all_intra_MOTU_dist_new <- unlist(intra_MOTU_dist_new)

#### Summary
cat('Mean intra-MOTU distance =', mean(all_intra_MOTU_dist_new), '\n',
    '95% CI = [', quantile(all_intra_MOTU_dist_new, 0.025), '-', quantile(all_intra_MOTU_dist_new, 0.975), ']\n',
    'Number of ASV pairs =', length(all_intra_MOTU_dist_new), '\n')


### 6-7-3. Check for barcoding gap
distance_df_new <- data.frame(distance = c(as.vector(inter_MOTU_dist_new),
                                           as.vector(all_intra_MOTU_dist_new)),
                              type = c(rep('Inter-MOTU', length(inter_MOTU_dist_new)),
                                       rep('Intra-MOTU', length(all_intra_MOTU_dist_new))))

### Plot histograms
barcoding_gap_after <- ggplot(distance_df_new, aes(x = distance, fill = type)) +
  geom_histogram(alpha = 0.5, position = 'identity', bins = 40) +
  scale_x_continuous(breaks = seq(0, max(distance_df$distance), by = 0.1)) +
  labs(x = 'Pairwise genetic distance (p-distance)',
       y = 'Frequency',
       fill = '') +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.3, 0.9),
        legend.justification = c(0.5, 1))
barcoding_gap_after

### Save in pdf
ggsave(barcoding_gap_after,
       file = './barcoding_gap/barcoding_gap_after.pdf',
       height = 10, width = 15, unit = 'cm')


## 6-8. Write MOTUs and ASVs tables
write.csv(ASV_table_updated, '../../data/ASV_table_updated.csv', row.names = F)    
write.csv(MOTU_table_updated, '../../data/MOTU_table_updated.csv', row.names = F)
# ___________________________________



# ---- 7. Occupancy modelling ----
## 7-1. Prepare data
### Import metadata
site_metadata <- read.csv2('../../data/site_metadata.csv')

### Create a long matrix with onr row per site per sample
site_long <- site_metadata %>%
  select(site = pool, spygen_code_1, spygen_code_2) %>%
  pivot_longer(cols = starts_with('spygen_code'),
               names_to = 'sample_col',
               values_to = 'sample_id') %>%
  group_by(site) %>%
  mutate(sample = factor(row_number(),
                         levels = c(1, 2),
                         labels = c('R1', 'R2'))) %>% 
  ungroup()

### Create a long matrix with one row per MOTU per sample
MOTU_long <- MOTU_table_updated %>%
  pivot_longer(cols = starts_with('SPY'),
               names_to  = 'sample_id',
               values_to = 'reads')

### Merge both long matrices
MOTU_site <- MOTU_long %>%
  left_join(site_long, by = c('sample_id')) %>%
  mutate(campaign = 'C1')

### Final 4D array
MOTU_array <- xtabs(reads ~ id + site + sample + campaign,
                    data = MOTU_site)
save(MOTU_array, file = './occupancy/MOTU_array.Rdata')


## 7-2. Run occupancy model
Y <- array(as.numeric(MOTU_array),
           dim = dim(MOTU_array),
           dimnames = dimnames(MOTU_array))
COI_occupancy <- Nemodel(array = Y,
                         covariates = NULL,
                         name = 'COI_occupancy',
                         protocol = 'seq_read',
                         nb_iterations = 50000,
                         nb_burnin = 25000,
                         nb_thinning = 5,
                         nb_chains = 3,
                         parallel = T)
save(COI_occupancy, file = './occupancy/COI_occupancy.Rdata')

### Check model convergence
model <- COI_occupancy@model
plot(model)

#### Gelman-Rubin statistic
Rhat <- model$BUGSoutput$summary[, 'Rhat']
hist(Rhat)
sum(Rhat > 1.1)
##### NOTE: The model converged successfully

#### Effective sample size
n_eff <- model$BUGSoutput$summary[, 'n.eff']
hist(n_eff)
sum(n_eff < 100)
#### NOTE: There is no autocorrelation between iterations


## 7-3. Calculate the minimum number of samples and minimum sequencing depth required to confidently detect MOTUs when present
COI_min_res <- min_resources(COI_occupancy,
                             resources = c('J', 'M'),
                             conf = 0.95)
save(COI_min_res, file = './occupancy/COI_min_res.Rdata')

### Link amplicons to scientific names
amplicons <- dimnames(MOTU_array)$id
amplicon_species <- MOTU_table_updated %>%
  select(id, ltg_name) %>%
  filter(id %in% amplicons) %>% 
  dplyr::rename(amplicon = id)

### 7-3-1. Filter MOTUs with insufficient J_min (= Minimum number of eDNA samples)
#### We build the data frame to record median and HDI of J_min for each MOTU
J_min_tab <- data.frame(amplicon = amplicons, # names of species
                        Median = COI_min_res@J_min$median[, 1, 1, 1],
                        Bottom = COI_min_res@J_min$hdi1[, 1, 1, 1],
                        Top = COI_min_res@J_min$hdi2[, 1, 1, 1]) %>%
  inner_join(amplicon_species)

#### We record MOTUs with J_min > 2 in the HDI
MOTU_to_remove_J <- filter(J_min_tab, Bottom > 2)
MOTU_to_remove_J
##### NOTE: 3 MOTUs should be removed here

### 7-3-2. Filter MOTUs with insufficient M_min (= Minimum sequencing depth)
#### Build the matrix with total number of reads per sample
M <- apply(MOTU_array, c(2, 3, 4), sum, na.rm = T)

#### Record MOTUs when both samples show M_min in a site
MOTU_to_remove_M <- data.frame(amplicon = character(),
                               Median = numeric(),
                               Bottom = numeric(),
                               Top = numeric(),
                               ltg_name = character())
for (n in dim(MOTU_array)[1]){
  for (i in dim(MOTU_array)[2]){
    comparison <- COI_min_res@M_min$hdi1[n, i, , ] >= M[i, , ]
    if(sum(comparison > 0)){
      row_to_add <- c(amplicon = amplicons[n],
                      Median = max(COI_min_res@M_min$median[n, i, , ]),
                      Bottom = max(COI_min_res@M_min$hdi1[n, i, , ]),
                      Top = max(COI_min_res@M_min$hdi2[n, i, , ])) %>% 
        inner_join(amplicon_species)
      MOTU_to_remove_M <- rbind(MOTU_to_remove_M,
                                row_to_add)        
    }
  }
}
MOTU_to_remove_M
##### NOTE: No MOTU to remove here

### 7-3-3. Remove the MOTUs for which we cannot ensure reliable detection
#### MOTU table
MOTU_table_reliable <- MOTU_table_updated %>% 
  filter(!id %in% MOTU_to_remove_J$amplicon)

#### ASV table
ASV_table_reliable <- ASV_table_updated %>% 
  filter(!MOTU %in% MOTU_to_remove_J$amplicon)


## 7-4. Write MOTUs and ASVs tables
write.csv(ASV_table_reliable, '../../data/ASV_table_reliable.csv', row.names = F)    
write.csv(MOTU_table_reliable, '../../data/MOTU_table_reliable.csv', row.names = F)
# ________________________________



# ---- 8. Accumulation curves ----
## 8-1. Create iNEXT input lists for abundance data
MOTU_array_reliable <- MOTU_array[!dimnames(MOTU_array)$id %in% MOTU_to_remove_J$amplicon, , , ]
sample_names <- c()
input_list <- list()
for (i in 1:dim(MOTU_array_reliable)[2]){
  for (j in 1:dim(MOTU_array_reliable)[3]){
    abund <- MOTU_array_reliable[, i, j]
    sample_names = c(sample_names,
                     paste0(colnames(MOTU_array_reliable)[i], '_sample', j))
    input_list <- c(input_list,
                    list(abund[abund > 0]))
  }
}
names(input_list) <- sample_names


## 8-2. Estimate accumulation curve for each sample
iNEXT_samples <- list()
for(i in 1:length(input_list)){
  iNEXT_samples <- append(iNEXT_samples,
                          list(iNEXT(input_list[i], q = 0, datatype = 'abundance', se = FALSE, endpoint = sum(input_list[[i]]))))
}
save(iNEXT_samples, file = './accumulations/iNEXT_samples.Rdata')


## 8-3. Combine results for all samples
iNEXT_all_samples <- map2_dfr(iNEXT_samples, names(input_list), function(obj, sample_name){
  df <- obj$iNextEst
  if (is.list(df) && !is.data.frame(df)){
    df <- df[[1]]
  } 
  df %>% mutate(Sample = sample_name)
})


## 8-4. Plot the curves
### 8-4-1. Plot sample curves independently
for(i in 1:length(input_list)){
  plot <- ggiNEXT(iNEXT_samples[[i]], type = 1) +
    labs(title = names(input_list)[i],
         x = 'Read count',
         y = 'MOTU richness') +
    theme_classic() +
    theme(legend.position = 'none')
  ggsave(plot,
         file = paste0('./accumulations/samples/', names(input_list)[i], '.pdf'),
         height = 7.5,
         width = 10,
         unit = 'cm')
}

### 8-4-2. Plot all curves
accumulation_samples_plot <- ggplot(iNEXT_all_samples,
                                    aes(x = m, # ou SC pour sampling coverage
                                        y = qD,
                                        color = Sample)) +
  geom_line(size = 1) +
  labs(x = 'Sampling coverage',
       y = 'MOTU richness') +
  theme_classic() +
  theme(legend.position = 'none')
accumulation_samples_plot

#### Save in pdf
ggsave(accumulation_samples_plot,
       file = './accumulations/accumulation_samples_plot_MOTU.pdf',
       height = 7.5,
       width = 12.5,
       unit = 'cm')


## 8-5. Check for sampling coverage
rarefaction_observed_values <- iNEXT_all_samples %>% 
  filter(Method == 'Observed')
print(rarefaction_observed_values$SC)
### NOTE: All samples reach a similar sampling coverage (> 0.999).
### There is no need to rarefy data.
# ________________________________



# ---- 9. Krona plot ----
### Correct invalid taxids
MOTUs_with_invalid_taxids <- MOTU_table_reliable %>% 
  filter(ltg_taxid < 0) %>% 
  select(ltg_name, ltg_taxid)
MOTUs_with_invalid_taxids

#### Manually check on NCBI taxonomy browser
####  - Oxydromus pallidus has taxid 3128520 (Oxydromus sp.)
####  - Spurilla neapolitana has taxid 929453
####  - Pavlova has taxid 2831
####  - Harpacticus gracilis has taxid 512228 (Harpacticus sp.)

MOTU_table_reliable$ltg_taxid[MOTU_table_reliable$ltg_name == 'Oxydromus pallidus'] <- 3128520
MOTU_table_reliable$ltg_taxid[MOTU_table_reliable$ltg_name == 'Spurilla neapolitana'] <- 929453
MOTU_table_reliable$ltg_taxid[MOTU_table_reliable$ltg_name == 'Pavlova'] <- 2831
MOTU_table_reliable$ltg_taxid[MOTU_table_reliable$ltg_name == 'Harpacticus gracilis'] <- 512228

### Create Krona input
Krona_tab <- data.frame()
for (i in 1:nrow(MOTU_table_reliable)){
  Krona_tab <- rbind(Krona_tab,
                     c(i, MOTU_table_reliable$ltg_taxid[i], 1))
}
colnames(Krona_tab) <- c('#query', '#taxID', '#score')
write_tsv(Krona_tab, file = './taxonomy/krona_input_mkLTG.tsv')

### Run Krona on a Linux terminal
'
\bash
conda activate krona
ktImportTaxonomy krona_input_mkLTG.tsv -o taxonomy.krona_mkLTG.html
'
# _______________________



# ---- 10. Non-indigenous MOTUs ----
## 10-1. Import Mediterranean NIS checklist from Galadini 2023 (https://www.mdpi.com/1424-2818/15/9/962)
nis_species_in_med <- read.csv2('./data/mediterranean_NIS.csv')


## 10-2. Add a 'Status' column in the MOTU table to record if the taxa is NIS
MOTU_table_final <- MOTU_table_reliable %>% 
  mutate(NIS = NA)
MOTU_table_final$NIS[MOTU_table_final$ltg_name %in% nis_species_in_med$Species] <- 'NIS'

### A literature screening revealed 3 other NIS:
###   - Amphibalanus amphitrite (Massé 2023, https://www.mdpi.com/1424-2818/15/2/161)
###   - Coryne eximia (Massé 2023, https://www.mdpi.com/1424-2818/15/2/161)
###   - Watersipora subtorquata (subatra) (Gauff 2023, https://aquaticinvasions.arphahub.com/article/108128/)

MOTU_table_final$NIS[MOTU_table_final$ltg_name == 'Amphibalanus amphitrite'] <- 'NIS'
MOTU_table_final$NIS[MOTU_table_final$ltg_name == 'Coryne eximia'] <- 'NIS'
MOTU_table_final$NIS[MOTU_table_final$ltg_name == 'Watersipora subtorquata'] <- 'NIS'

### Schizoporella pungens is a potential NIS ! Not recorded yet in the Mediterranean waters
### but established in Madeira Island (Canning-Clode 2013, http://www.aquaticinvasions.net/2013/issue3.html)

MOTU_table_final$NIS[MOTU_table_final$ltg_name == 'Schizoporella pungens'] <- 'NIS'

## 10-3. Write MOTUs and ASVs tables
ASV_table_final <- ASV_table_reliable
write.csv(ASV_table_final, '../../data/ASV_table_final.csv', row.names = F)    
write.csv(MOTU_table_final, '../../data/MOTU_table_final.csv', row.names = F)
# __________________________________