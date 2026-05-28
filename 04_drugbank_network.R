#DrugBank Analysis
# CRAN
install.packages("remotes")
remotes::install_github("ropensci/dbparser")
library(dbparser)


library(dbparser)
library(dplyr)
library(readr)
library(xml2)


#parse the drugbank

   dv = parseDrugBank(
    db_path =  ("/Users/ebraraltinalan/Desktop/GEO/full_database.xml"),
     drug_options = drug_node_options(),
     references_options = references_node_options(),
     cett_options = cett_nodes_options(),  # carriers, enzymes, targets, transporters
     parse_products = TRUE,
     parse_salts = TRUE
   )
   
   str(dv, max.level = 1) #check parseDrugBank()
   
   #List of 5
   #$ drugs     :List of 28
   #$ salts     : tibble [2,922 × 8] (S3: tbl_df/tbl/data.frame)
  # $ products  : tibble [455,970 × 19] (S3: tbl_df/tbl/data.frame)
   #$ references:List of 5
   #$ cett      :List of 4
   #- attr(*, "class")= chr "dvobject"
   #- attr(*, "original_db_info")=List of 3
   
   
   names(dv$cett)
   names(dv$cett$targets)
   str(dv$cett$targets$polypeptides, max.level = 2)
   
   
   # Extract target and polypeptide data
   targets <- dv$cett$targets
   polypeptides <- targets$polypeptides$general_information
   synonyms <- dv$cett$targets$polypeptides$synonyms
   
   # Replace gene symbols
   symbol_map <- setNames(Homo_sapiens_gene_info$V3, Homo_sapiens_gene_info$V5)
   polypeptides$gene_name <- ifelse(
     polypeptides$gene_name %in% names(symbol_map),
     symbol_map[polypeptides$gene_name],
     polypeptides$gene_name
   )
   
   targets_general <- targets$general_information  # tibble
   nrow(targets_general)
   
   drug_names <- drugs_general[, c("drugbank_id", "name")]
   
   
   nrow(drugs_general)
   nrow(synonyms)
   nrow(targets)
   nrow(polypeptides)
   
   #merge tables 
   library(dplyr)
   
   # The targets will merge with polypeptides (by target_id)
   targets_merged <- targets_general %>%
     left_join(polypeptides, by = "target_id") 
   
   # Merge targets + polypeptides with synonym (by target_id)
   targets_merged <- targets_merged %>%
     left_join(synonyms, by = "target_id")     
   
   
   # Add drugs informations (by drugbank_id)
   final_merged <- targets_merged %>%
     left_join(drugs_general[, c("drugbank_id", "name")], by = "drugbank_id")

   # Filter   the “Known action”  interactions
   final_verified <- final_merged %>%
     filter(known_action == "yes")

   
   # Control
   nrow(final_verified)
   head(final_verified$known_action)
   
   library(readxl)

   

   meta <- read_excel("/Users/ebraraltinalan/Desktop/GEO/final-nocov--GEO/analysis/meta_analysis_sig_genes_names.xlsx")
   
   # Control
   head(meta)
   colnames(meta)
   
   #categorize by up-regulated and down regulated
   meta <- meta %>%
     mutate(regulation = case_when(
       meta_LFc > 0  ~ "up-regulated",
       meta_LFc < 0  ~ "down-regulated",
       TRUE          ~ "no_change"
     ))
   
   
   colnames(final_verified)

   
   
   
   
   
   library(dplyr)
   
   # ️Change gene names to uppercase
   meta <- meta %>%
     mutate(`Gene Symbol` = toupper(`Gene Symbol`)) 
   
   final_verified <- final_verified %>%
     mutate(gene_name = toupper(gene_name))
   
   final_verified_unique <- final_verified %>%
     distinct(gene_name, drugbank_id, .keep_all = TRUE)
   
   
   # Merge: only matching genes
   drug_gene_interactions <- inner_join(
     meta,
     final_verified_unique,
     by = c("Gene Symbol" = "gene_name")
   )
   head(drug_gene_interactions)
   
   library(writexl)
   write_xlsx(drug_gene_interactions,"drug_gene_interactions_unique.xlsx")
   
   
   drug_gene_interactions_unique <- drug_gene_interactions %>%
     select(
       `Gene Symbol`,
       GeneName,
       regulation,
       meta_LFc,
       meta_pval,
       drugbank_id,
       name,           # Drug primary name
       target_id,
       known_action,
       general_function,
       specific_function,
       synonym
     )
   
   write_xlsx(drug_gene_interactions_unique, "drug_gene_interactions_unique_filtered.xlsx")
   
   # Control
   head(drug_gene_interactions)
   nrow(drug_gene_interactions)
   
   
   #Visualisation
  
   #build drug-gene network (igraph + ggraph)
 
   library(igraph)
   library(ggraph)
   library(ggplot2)  
   
   # Create edge list
   edges = drug_gene_interactions_unique_filtered[, c('Gene Symbol', 'Drug name', "general_function")]  # gene -> drug
   colnames(edges) = c("from", "to", "general_function")

   genes = unique(drug_gene_interactions_unique_filtered[, c('Gene Symbol', "meta_LFc")])
   genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
   colnames(genes)[1] = "name"  # for igraph
   
   # Create a unique list of drug nodes
   drugs = unique(data.frame(name = drug_gene_interactions_unique_filtered$`Drug name`, regulation = "drug"))
   drugs$meta_LFc <- NA
   
  
   # Combine all nodes
   nodes = rbind(genes, drugs)  
   #build i graph
   graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
   
   # Plot with ggraph
   ggraph(graph, layout = "fr") +
     geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
     geom_node_point(aes(color = regulation), size = 14) +  #  nodes
     geom_node_text(aes(label = name), size = 5, color = "black") +  # label centered
     scale_color_manual(
       values = c(
         "up" = "#FF9999",      
         "down" = "#9999FF",
         "drug" = "#A8D5BA"
       )
     ) +
     theme_void()
   
   
 
   library(dplyr)
   
   top_genes <- drug_gene_interactions_unique_filtered %>%
     arrange(desc(abs(meta_LFc))) %>%
     distinct(`Gene Symbol`, .keep_all = TRUE) %>%
     slice(1:15)
   
   # top15
   edges_top <- drug_gene_interactions_unique_filtered %>%
     filter(`Gene Symbol` %in% top_genes$`Gene Symbol`) %>%
     select(`Gene Symbol`, `Drug name`, general_function)
   
   colnames(edges_top) <- c("from", "to", "general_function")
   
   # gene nodes
   genes_top <- unique(top_genes[, c('Gene Symbol', "meta_LFc")])
   genes_top$regulation <- ifelse(genes_top$meta_LFc > 0, "up", "down")
   colnames(genes_top)[1] <- "name"
   
   # drug nodes
   drugs_top <- unique(data.frame(
     name = edges_top$to,
     regulation = "drug",
     meta_LFc = NA
   ))
   
   # bind nodes
   nodes_top <- rbind(genes_top, drugs_top)
   
   # Graphs
   graph_top <- graph_from_data_frame(d = edges_top, vertices = nodes_top, directed = TRUE)
   
   # 
   ggraph(graph_top, layout = "fr") +
     geom_edge_link(alpha = 0.8, color = "gray70") +
     geom_node_point(aes(color = regulation), size = 8) +
     geom_node_text(aes(label = name), size = 4, color = "black") +
     scale_color_manual(values = c("up" = "#FF9999", "down" = "#9999FF", "drug" = "#A8D5BA")) +
     theme_void()
   

   
#Various graphic designs were tested in separate scripts
   
   drug_names = dv$drugs[[1]][, c("drugbank_id", "name")]
   colnames(drug_names)[2] = "drug_name"

   
   #classify the drugs
   install.packages(c("webchem", "dplyr", "purrr", "readr"))
   install.packages("PubChemR")
   library(PubChemR)
   library(webchem)
   library(dplyr)
   library(purrr)
   
   #convert drug names into PubChem CIDs
   drug_names_dge = drug_gene_interactions_unique_filtered$`Drug name`
   cids = get_cid(drug_names_dge, from = "name", match = "first", domain = "compound")
   
   #remove rows with NA
   cids = cids[!is.na(cids$cid),]
   cids = cids[!grepl("^4", cids$query), ]
   cids = cids[!grepl("^6", cids$query), ]
   cids = cids[!grepl("^3", cids$query), ]
   
   write(cids$cid, "cids.txt")
   
   category = dv$drugs$categories
   
   drugs_with_categories = merge(
     drug_names,
     category,
     by.x = "drugbank_id",
     by.y = "drugbank_id"
   )
   
   drugs_useful = drugs_with_categories[drugs_with_categories$drug_name %in% cids$query,]
   
   #remove rows with missing mesh_id
   drugs_useful$mesh_id[drugs_useful$mesh_id == ""] = NA
   drugs_useful = drugs_useful[!is.na(drugs_useful$mesh_id),]
   
   #remove rows with non-pharmacological/medical classification
   not_category = c("Hormones", "Proteins", "Peptides", 
                    "Enzyme Inhibitors", "Benzazepines", "Benzodiazepinones", "Nervous System",
                    "Peripheral Nervous System Agents", "Alcohols", "Central Nervous System Agents",
                    "Central Nervous System Depressants", "Gastrointestinal Agents", "Chlorohydrins",
                    "Amines", "Autonomic Agents", "Neurotransmitter Agents",
                    "Peripheral Nervous System Agents", "Phenethylamines", "Cardiovascular Agents",
                    "Hematologic Agents", "Purinergic Agents", "Pyridines", "Sulfur Compounds",
                    "Thiophenes", "Ethers", "Methyl Ethers", "Amines", "Oligopeptides",
                    "Pyrimidines", "Pyrimidinones", "Amides", "Hexoses", "Ketoses",
                    "Imidazoles", "Acids, Acyclic", "Acids, Carbocyclic",
                    "Acetates", "Adenine Nucleotides", "Adrenal Cortex Hormones",
                    "Adrenergic Agents", "Amino Acids", "Amino Acids, Peptides, and Proteins",
                    "Anesthetics, Inhalation", "Anesthetics, Intravenous", "Anesthetics, General",
                    "Anemia, Iron-Deficiency", "Analgesics, Non-Narcotic", "Aminobutyrates",
                    "Alkaloids", "Alkenes", "Aniline Compounds", "Benzene Derivatives", "Benzodiazepines and benzodiazepine derivatives",
                    "Benzopyrans", "Biological Factors", "Carbohydrates", "Coenzymes", "Cytochrome P-450 CYP1A2 Inducers",
                    "Heterocyclic Compounds, Fused-Ring", "Cytochrome P-450 CYP1A2 Inhibitors", 
                    "Cytochrome P-450 CYP2B6 Inhibitors", "Cytochrome P-450 CYP2C19 Inhibitors",
                    "Cytochrome P-450 CYP2C8 Inhibitors", "Cytochrome P-450 CYP2C9 Inhibitors",
                    "Cytochrome P-450 CYP2D6 Inhibitors", "Cytochrome P-450 CYP2E1 Inhibitors",
                    "Cytochrome P-450 Enzyme Inhibitors", "GABA Agents", "GABA Modulators",
                    "Cytochrome P-450 CYP3A Inducers", "Cytochrome P-450 CYP3A4 Inducers",
                    "Cytochrome P-450 Enzyme Inducers", "Monosaccharides", "Adrenergic alpha-2 Receptor Agonists",
                    "Hydroxymethylglutaryl-CoA Reductase Inhibitors", "Alkylating Drugs", "Carbamates",
                    "Eicosanoids", "Lipids", "Fatty Acids", "Vasodilation", "Purinergic P2 Receptor Antagonists",
                    "Purinergic Antagonists", "Nitroso Compounds", "Nitrosourea Compounds", "GABA-A Receptor Agonists",
                    "Muscle Relaxants, Centrally Acting Agents", "Fused-Ring Compounds", "Digoxin and derivatives",
                    "Indoles", "Indole Alkaloids", "Purines", "Ribonucleotides", "Purine Nucleotides", "Nucleic Acids, Nucleotides, and Nucleosides",
                    "Protein Kinase Inhibitors", "Tyrosine Kinase Inhibitors", "Iron Compounds",
                    "Depression, Postpartum", "Organometallic compounds", "Minerals", "Neuromuscular Agents",
                    "Delayed-Action Preparations", "Physiological Phenomena", "Diet, Food, and Nutrition", "Pyrazoles",
                    "Octanols", "Nicotinic Acids", "Pyrrolidines", "Tumor Suppressor Proteins", "Starch", 
                    "Corpus Luteum Hormones", "Cyclodextrins", "Dextrins", "Dietary Carbohydrates", 
                    "Glucans", "Polysaccharides", "Pregnenediones", "Pregnenes", "Progesterone Congeners",
                    "Anticholesteremic Agents", "Noxae", "Photosensitizing Agents", "Toxic Actions" , 
                    "Gastrins", "Gastrointestinal Hormones", "Hormones, Hormone Substitutes, and Hormone Antagonists",         
                    "Serotonergic Drugs Shown to Increase Risk of Serotonin Syndrome", "Purinergic P2Y Receptor Antagonists",                            
                    "Thienopyridines", "Nicotinic Antagonists", "Antineoplastic Agents", "Antineoplastic Agents, Alkylating",                              
                    "Adrenergic Agonists", "Adrenergic alpha-Agonists", "Adjuvants, Anesthesia",                                          
                    "GABA Agonists", "Antihypertensive Agents", "Autacoids", "Fatty Acids, Unsaturated",                                       
                    "Prostaglandins I", "Compounds used in a research, industrial, or household setting", 
                    "Digitalis Glycosides","Protective Agents","Piperazines", "Pyrazines", "Biogenic Amines",                                                
                    "Biogenic Monoamines", "Catechols","Ethanolamines","Cytochrome P-450 CYP2C19 Inducers", 
                    "Cytochrome P-450 CYP2C8 Inducers",                               
                    "Cytochrome P-450 CYP2C9 Inducers", "Sleep Aids, Pharmaceutical", "Cytochrome P-450 CYP2B6 Inducers",                               
                    "Cytochrome P-450 CYP3A Inhibitors", "Cytochrome P-450 CYP3A4 Inhibitors", "Antidotes",
                    "Biological Products", "Central Nervous System Stimulants", "Complex Mixtures", "Cyclohexanes",                                                  
                    "Cycloparaffins", "GABA Antagonists", "GABA-A Receptor Antagonists", "Lactones", "Pharmaceutical Preparations",                                    
                    "Plant Extracts", "Plant Preparations", "Sesquiterpenes", "Terpenes", "Toxins, Biological",                                             
                    "Cytosine Nucleotides", "Nucleotides", "Pyrimidine Nucleotides", "Dermatologicals",                                                
                    "Pregnadienes", "Pregnanes", "Steroids, Fluorinated", "Thiobarbiturates",                                                
                    "Naphthalenes", "Hydrocarbons, Halogenated", "Cytochrome P-450 CYP2E1 Inducers",                               
                    "Phenobarbital and similars", "Phenols",                                                        
                    "Cystine Depleting Agents", "Ethylamines", "Mercaptoethylamines","Sulfhydryl Compounds",                                           
                    "Anti-Infective Agents, Local", "Solvents","Cyclooxygenase Inhibitors", "Hydroxybenzoates",                                               
                    "Salicylates",                                                    
                    "Sensory System Agents", "Ethyl Ethers","Heptanoic Acids", "Cytochrome P-450 CYP3A5 Inhibitors",                             
                    "Hydrocarbons, Fluorinated", "Sulfones","Histamine Agents", "Histamine H1 Antagonists",                                       
                    "Histamine H1 Antagonists, Non-Sedating", "Piperidines", "Chemically-Induced Disorders",                                   
                    "Cytochrome P-450 CYP3A5 Inducers",                               
                    "Cytochrome P-450 CYP3A7 Inducers",                               
                    "Excitatory Amino Acid Agents",                                   
                    "Excitatory Amino Acid Antagonists",                             
                    "Aza Compounds",                                                 
                    "Sleep Initiation and Maintenance Disorders",                     
                    "Piperidones",                                                    
                    "Alkanes",                                                        
                    "Alkanesulfonic Acids",                                           
                    "Hydrocarbons, Acyclic",                                          
                    "Sulfonic Acids",                                                 
                    "Sulfur Acids",                                                   
                    "Enzymes and Coenzymes" ,                                         
                    "Prostaglandins D",                                               
                    "Flavins"    ,                                                    
                    "Pigments, Biological"  ,                                         
                    "Pteridines"         ,                                            
                    "Antimetabolites"        ,                                        
                    "Dicarboxylic Acids"     ,                                        
                    "Glutarates"             ,                                        
                    "Fatty Acids, Volatile"  ,                                        
                    "Valerates"              ,                                        
                    "17-Ketosteroids"         ,                                       
                    "Gonadal Steroid Hormones"  ,                                     
                    "Bridged-Ring Compounds" ,                                        
                    "Flavones" ,                                                      
                    "Flavonoids",                                                     
                    "Pyrans" ,                                                        
                    "Nucleosides",                                                    
                    "Purine Nucleosides" ,                                            
                    "Ribonucleosides"    ,                                            
                    "Aminopyridines"   ,                                              
                    "Hydroxy Acids"    ,                                              
                    "Immunologic Factors" ,                                           
                    "Dioxoles"       ,                                                
                    "Phosphodiesterase Inhibitors" ,                                  
                    "Ammonium Compounds"   ,                                          
                    "Nitrogen Compounds"   ,                                          
                    "Onium Compounds"       ,                                         
                    "Surface-Active Agents"  ,                                        
                    "Abortifacient Agents"   ,                                        
                    "Hormonal Contraceptives for Systemic Use"      ,                 
                    "Luteolytic Agents"   ,                                           
                    "Prostaglandins F, Synthetic"     ,                               
                    "Prostaglandins, Synthetic"   ,                                   
                    "Reproductive Control Agents"   ,                                 
                    "Oxazines",                                         
                    "Cytochrome P-450 CYP3A7 Inhibitors"   ,                          
                    "Carbon Radioisotopes",                                 
                    "Fatty Alcohols" ,                                    
                    "Hexanols",                                        
                    "Food",                                        
                    "Micronutrients",                                                                                            
                    "Abortifacient Agents, Nonsteroidal",                          
                    "Prostaglandins F",                                               
                    "Uterotonic agents", "Organometallic Compounds","Aminobenzoates", "Thiazoles", "Hydrazines" )
   drugs_useful = drugs_useful[!drugs_useful$category %in% not_category,]
   
   print(unique(drugs_useful$category))
   unique(drugs_useful$drug_name)
   
   
   drugs_useful_grouped = drugs_useful %>%
     group_by(drug_name) %>%
     summarise(
       drugbank_id = first(drugbank_id),  
       category = list(unique(category)), # collect unique categories
       mesh_id = list(unique(mesh_id))    # collect unique MeSH IDs
     )
   drugs_useful_grouped = drugs_useful_grouped[, 1:3]
   
   print(unique(drug_gene_interactions_unique_filtered$`Drug name`))
   
   devtools::install_github("yduan004/drugbankR")
   library(drugbankR)
   
   drugs = unique(drug_gene_interactions_unique_filtered$`Drug name`)

  