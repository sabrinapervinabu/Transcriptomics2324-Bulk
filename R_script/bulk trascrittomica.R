#### DOWNLOAD AND INSTALL PACKAGES AND LIBRARIES ####
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("recount3")

if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
  "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
  "derfinder", "GenomicState", "bumphunter", "derfinderPlot", "sessioninfo"
))

#da runnare separatamente da tutto il codice di installazione 
library(recount3)
library(recount)

library(edgeR)

####INCUBO####
#rse_brain <- readRDS("rse_brain.RDS")
#assays(rse_brain)$counts <- transform_counts(rse_brain)
#assays(rse_brain)$TPM <- recount::getTPM(rse_brain)

#### UPLOAD THE TISSUES WITH ALL THE SAMPLES AND THEIR COUNTS of every GENE ####
rse_brain <- readRDS("rse_brain.RDS")
rse_heart <- readRDS("rse_heart.RDS")
rse_kidney <- readRDS("rse_kidney.RDS")
# the counts are in “coverage” format, not count of reads

# To transform them in counts as we need for DE analysis:
assays(rse_brain)$counts <- transform_counts(rse_brain)
rse_brain
assays(rse_heart)$counts <- transform_counts(rse_heart)
rse_heart
assays(rse_kidney)$counts <- transform_counts(rse_kidney)
rse_kidney

#### SELEZIONE DELLE COLONNE PER RIN, %rRNA & %MAPPED READS - 3 SAMPLES PER TESSUTO ####
# What is the minimum RIN? I usually recommend at least 7, 
# but 6 or higher is usually considered to be “acceptable”. So all the three samples pass the threshold.
colData(rse_brain)$gtex.smrin[30] #poi lo fa per le altre colonne anche
colData(rse_brain)$gtex.smrin[38]
colData(rse_brain)$gtex.smrin[41]

# Estimated fraction of rRNA: #This should be very low, never anyway higher than 10% 
# (or 0.1 since here is the fraction to be reported). Once again, all three samples pass the threshold.
colData(rse_brain)$gtex.smrrnart[30] #poi lo fa per le altre colonne anche
colData(rse_brain)$gtex.smrrnart[38]
colData(rse_brain)$gtex.smrrnart[41]

# Finally, the percentage of mapped reads: 
# we want here at least 85% of the reads uniquely mapped, 
# since it is a human sample
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[30] #poi lo fa per le altre colonne anche
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[38]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[41]

rse_brain_selected <- rse_brain[,c(30,38,41)]
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_brain_selected

colData(rse_heart)$gtex.smrin[30]
colData(rse_heart)$gtex.smrin[32]
colData(rse_heart)$gtex.smrin[33]

colData(rse_heart)$gtex.smrrnart[30] 
colData(rse_heart)$gtex.smrrnart[32]
colData(rse_heart)$gtex.smrrnart[33]

colData(rse_heart)$"recount_qc.star.uniquely_mapped_reads_%_both"[30] 
colData(rse_heart)$"recount_qc.star.uniquely_mapped_reads_%_both"[32] 
colData(rse_heart)$"recount_qc.star.uniquely_mapped_reads_%_both"[33] 

rse_heart_selected <- rse_heart[,c(30,32,33)]
counts_heart_selected <- assays(rse_heart_selected)$counts
counts_heart_selected 

colData(rse_kidney)$gtex.smrin[49]
colData(rse_kidney)$gtex.smrin[56]
colData(rse_kidney)$gtex.smrin[63]

colData(rse_kidney)$gtex.smrrnart[49] 
colData(rse_kidney)$gtex.smrrnart[56]
colData(rse_kidney)$gtex.smrrnart[63]

colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[49] 
colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[56] 
colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[63] 

rse_kidney_selected <- rse_kidney[,c(49,56,63)]
counts_kidney_selected <- assays(rse_kidney_selected)$counts
counts_kidney_selected 

#### FORMAZIONE DELLA TABELLA FINALE (9 SAMPLES) + FORMIAMO OBJECT DGE x edgeR ####
# We build the count table, give more meaningful name to the columns and the rows, 
# and finally build the DGE object for edgeR.
x <- cbind(counts_brain_selected,counts_heart_selected, counts_kidney_selected)
colnames(x) <- c("Brain30", "Brain38","Brain41","Heart30", "Heart32","Heart33",
                 "Kidney49","Kidney56","Kidney63")
rownames(x) <- rowData(rse_brain_selected)$gene_name #diamo il nome dei geni ad ogni riga 

y <- DGEList(counts=x)

# We define how replicates are grouped:
group <- as.factor(c("Brain","Brain","Brain","Heart","Heart","Heart","Kidney","Kidney","Kidney"))
y$samples$group <- group 
    #così edgeR capisce quali a quale gruppo appartiene ogni sample - e confronta i gruppi

# We add to the samples the “quality” information that we employed to select them. 
# For this we can add new labels to the “samples” field of the DGE object, 
# like we just did for the groups.
y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,
                             colData(rse_heart_selected)$gtex.smrin,
                             colData(rse_kidney_selected)$gtex.smrin)) 
                            #ad ogni sample si assegna il RIN

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,
                               colData(rse_heart_selected)$gtex.smtsd,
                               colData(rse_kidney_selected)$gtex.smtsd))
                            #ad ogni sample si assegna la sezione del tessuto da cui è stato prelevato 

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,
                             colData(rse_heart_selected)$gtex.sex,
                             colData(rse_kidney_selected)$gtex.sex))
                            #ad ogni sample si assegna il sesso del soggetto 

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,
                             colData(rse_heart_selected)$gtex.age,
                             colData(rse_kidney_selected)$gtex.age))
                            #ad ogni sample si assegna l'età

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,
                              colData(rse_heart_selected)$gtex.smrrnart,
                              colData(rse_kidney_selected)$gtex.smrrnart))
                            #ad ogni sample si assegna la %rRNA 

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", 
                                colData(rse_heart_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",
                                colData(rse_kidney_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))
                            

#ad ogni sample si assegna la %uniquely mapped reads 

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", 
                              colData(rse_heart_selected)$"recount_qc.aligned_reads%.chrm",
                              colData(rse_kidney_selected)$"recount_qc.aligned_reads%.chrm"))
                            #ad ogni sample si assegna la %reads mapped on Cromosoma Mitocondriale 

y
dim(y) 

#The combined count data is organized into a DGEList object, y, which is the standard input format for 
#the edgeR package. This object includes the counts matrix and metadata such as group (tissue type), 
#RIN, tissue section, sex, age, rRNA percentage, uniquely mapped reads percentage, 
#and mitochondrial chromosome reads percentage. 
#This metadata is essential for understanding the potential sources of variability in the data.


# Il risultato è una tabella che mostra il numero di geni che hanno zero conteggi 
# in tutti i campioni (TRUE) e il numero di geni che non hanno zero conteggi in tutti i campioni (FALSE).
table(rowSums(y$counts==0)==9) #righe, e quindi geni, dove tutte e 9 i campioni hanno 0 reads 
        #FALSE  TRUE 
        #39510 14532 
# A lot of genes with zero counts (i TRUE).

keep.exprs <- filterByExpr(y, group=group)
  # La funzione filterByExpr viene utilizzata per identificare i geni che hanno un'espressione sufficiente 
  # in almeno uno dei gruppi definiti nel vettore group. 
  # Questo passaggio è importante per rimuovere i geni con espressione molto bassa o nulla, 
  # che possono introdurre rumore nell'analisi.

y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  # Questa linea applica il filtraggio ai dati di conteggio nell'oggetto y, 
  # mantenendo solo i geni che hanno passato il filtro (keep.exprs). 
  # L'opzione keep.lib.sizes=FALSE indica che le dimensioni delle librerie 
  # (totale dei conteggi per campione) non devono essere mantenute.

# A lot of genes removed because of zero or low expression. 

dim(y) #restituisce le dimensioni del dataset filtrato, 
      # mostrando il numero di geni rimanenti e il numero di campioni.

# edgeR includes several useful functions for transforming the counts into counts per million, 
# FPKM, TPM etc. Let us extract and store in a vector the log of the counts per million 
# before normalization with the “cpm” function, and then normalize them:
logcpm_before <- cpm(y, log=TRUE)
      # La funzione cpm calcola i counts per million (CPM) e, con l'opzione log=TRUE, 
      # restituisce i logaritmi dei CPM. Questi valori vengono salvati in logcpm_before.
      # Salvare i log CPM prima della normalizzazione permette di confrontare 
      # l'effetto della normalizzazione sui dati.

# La funzione calcNormFactors calcola i fattori di normalizzazione per i dati di conteggio 
# utilizzando il metodo TMM (Trimmed Mean of M-values). 
y <- calcNormFactors(y, method = "TMM")
y 

    # I "fattori di normalizzazione" si riferiscono ai coefficienti o ai parametri calcolati e 
    # applicati ai dati di espressione genica per compensare le variazioni tecniche tra campioni 
    # nei dati di RNA-Seq. Questi fattori sono essenziali per rendere i dati comparabili tra campioni 
    # e per garantire che le analisi successive siano accurate e robuste

# Notice that all of the kidney samples have a very small normalization factor. 
# Let us compare the distribution of count values before and after normalization:
logcpm_after <- cpm(y, log=TRUE)
boxplot(logcpm_before, notch=T, las = 2, main = "before TMM") 

boxplot(logcpm_after, notch=T, las = 2, main = "after TMM")
    # post normalizzazione, sono stati corretti gli errori sistematici - per cui si finisce per 
    # medie di questi CPM più uniformi e maggiormente comparabili tra loro 

# Now we design the linear model. Intercept or not, it changes little. 
# Anyway, from a logical point of view the intercept is not needed here.
  # model.matrix() è una funzione in R che crea una matrice di progettazione per modelli lineari multipli, 
        # partendo da una specifica formula.
  # ~0+group specifica la formula del modello lineare. Il 0+ indica di non includere 
        # un termine di intercetta nel modello, quindi il modello sarà specificamente 
        # basato sui gruppi definiti nella variabile group.
  # data=y$samples specifica che i dati per costruire il modello sono 
        #contenuti nel campo samples dell'oggetto y
design <- model.matrix(~0+group, data=y$samples) #crea un modello lineare basato sui gruppi 
colnames(design) <- levels(y$samples$group) # assegna nome alle colonne 
design                                      # Questo comando visualizza la matrice di progettazione


# Let’s now see if/how the samples/replicates cluster together.
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group, main= 'Tissue')
#Replicates cluster well. BUT...

# ... One of the brain samples is a bit distant from the other two. 
# By trying labelling the points in the plot with different info, 
# I could single out what seem to be the most relevant factors:
      # si eseguono grafici MDS (Multi-Dimensional Scaling) 
      # con l'obiettivo è visualizzare come i campioni (o replicati) 
      # si raggruppano o differenziano in base alla loro espressione genica
plotMDS(logcpm, labels=y$samples$rRNA, main= '%rRNA' ) 

plotMDS(logcpm, labels=y$samples$chrm, main= 'chrM')

plotMDS(logcpm, labels=y$samples$age, main= 'age')

# That is - the more distant brain replicate differs from the other two 
# for the % of rRNA and of mitochondrial RNA -- and NOT the age of the donor.


y <- estimateDisp(y, design) # La dispersione è una misura della variabilità nell'espressione genica 
                          # tra i campioni, che può includere variabilità biologica e tecnica.
      # design è la matrice di progettazione che specifica il modello lineare utilizzato per l'analisi 

plotBCV(y, main= 'BVC') #visualizza il grafico del coefficiente di variazione biologica (BCV) stimato per ogni gene. 
  # Il BCV è una misura della variabilità relativa tra i campioni, 
  # calcolata come il rapporto tra la deviazione standard e il valore medio per un gene.

  # Un BCV alto indica una maggiore variabilità tra i campioni per quel gene.

# La linea "Common" nel grafico plotBCV rappresenta il BCV comune stimato per tutti i geni nel dataset.
  # Un singolo valore stimato che indica la variabilità media tra i geni nel dataset. 
  # Indica quanto la varianza tra i geni è comune a tutti, 
  # fornendo un riferimento per valutare la variabilità relativa di ciascun gene.

# La linea "Trend" nel grafico plotBCV rappresenta il trend di variazione del BCV stimato VS conteggio medio.
  # Questo trend indica se il BCV stimato tende a variare in relazione ai livelli di espressione genica. 
  # Ad esempio, potrebbe mostrare se il BCV è più alto o più basso per i geni con bassi 
  # livelli di espressione rispetto a quelli con alti livelli di espressione.

# Ogni punto nel grafico rappresenta un gene nel dataset. L'asse x di solito rappresenta il conteggio medio 
# o altre misure di espressione, mentre l'asse y rappresenta il BCV stimato per quel gene.

#I punti che si trovano sopra la linea "Common" hanno un BCV superiore alla media comune stimata, 
# indicando una maggiore variabilità tra i geni rispetto alla media.
# Il trend generale dei punti rispetto alla linea "Trend" fornisce un'indicazione se il BCV tende 
# a cambiare in modo sistemico con l'espressione genica.

# Quite high variability (what is the BCV on the y axis, 
# and how does it relate to the NB distribution employed to model the data?)
fit <- glmQLFit(y, design)
    # glmQLFit() stimola il modello di regressione log-lineare utilizzando i dati di dispersione stimati (y) 
    # e la matrice di progettazione (design)

#heart (top) vs brain (bottom)
qlfHB <- glmQLFTest(fit, contrast=c(-1,1,0)) #Definizione del test di contrasto specifico.
qlfBH  <- glmQLFTest(fit, contrast=c(1,-1,0))
  # The three comparisons return a data structure 
  # with all the info regarding the comparisons themselves, 
  # including a table with the result of the test for each gene.

  #“topTags” extracts the table, sorted by p-value, and adds the FDR:
topTags(qlfBH, n=10,adjust.method = "BH", sort.by = "PValue") 
  #Estrazione dei geni con significative differenze di espressione, ordinati per p-value corretto.

#logFC (Log Fold Change): Rappresenta la differenza media nei logaritmi dei conteggi 
#tra i campioni di heart e Brain. 
# logFC > 0 il gene è sovraespresso nei campioni di heart rispetto a quelli di Brain
# logFC < 0 indica sottoespressione nei campioni di heart rispetto a quelli di Brain.

#logCPM (Log Counts Per Million): Media del logaritmo dei conteggi per milione di reads (CPM) per il gene.

#F (Statistiche del Test): Statistica del test (F-statistic) 
 #che valuta la significatività della differenza di espressione tra i due gruppi.

# FDR (False Discovery Rate): Correzione del valore di p per il controllo del tasso di falsi positivi 
# utilizzando il metodo di correzione multiplo di Benjamini-Hochberg (BH). 
# Un valore di FDR più basso indica una maggiore significatività statistica dei risultati.


# CNTNAP4: Ha un logFC positivo di 9.818029, indicando un'elevata sovraespressione nei campioni di BRAIN 
# rispetto a quelli di HEART. Il valore di p (PValue) è molto basso (3.217552e-08), 
# con un FDR altrettanto basso (0.0002803556), suggerendo che questa differenza di espressione 
# è altamente significativa statisticamente.

###### RIPRENDI DA QUI #####
qlfKB <- glmQLFTest(fit, contrast=c(-1,0,1))
qlfKH <- glmQLFTest(fit, contrast=c(0,-1,1))
qlfHB <- glmQLFTest(fit, contrast=c(-1,1,0))
qlfHK <- glmQLFTest(fit, contrast=c(0,1,-1))
qlfBH <- glmQLFTest(fit, contrast=c(-1,-1,0))
qlfBK <- glmQLFTest(fit, contrast=c(-1,0,-1))
summary(decideTests(qlfKB, p.value=0.01, adjust.method = "BH",lfc=1))
summary(decideTests(qlfKH, p.value=0.01, adjust.method = "BH",lfc=1))
summary(decideTests(qlfHB, p.value=0.01, adjust.method = "BH",lfc=1))
summary(decideTests(qlfHK, p.value=0.01, adjust.method = "BH",lfc=1))
summary(decideTests(qlfBH, p.value=0.01, adjust.method = "BH",lfc=1))
summary(decideTests(qlfBK, p.value=0.01, adjust.method = "BH",lfc=1))


#kidney (top) vs brain (bottom)
qlfKB <- glmQLFTest(fit, contrast=c(-1,0,1))
topTags(qlfKB, n=10,adjust.method = "BH", sort.by = "PValue")
# ITGB6: Ha un logFC positivo di 12.460020, indicando un'elevata sovraespressione nei campioni di kidney 
# rispetto a quelli di Brain. Il valore di p (PValue) è molto basso (9.943291e-09), 
# con un FDR altrettanto basso (0.0002296701), suggerendo che questa differenza di espressione è altamente 
# significativa statisticamente

#kidney (bottom) vs heart (top)
qlfHK <- glmQLFTest(fit, contrast=c(0,1,-1))
topTags(qlfHK, n=10,adjust.method = "BH", sort.by = "PValue")
# KCNJ16 ha un logFC NEGATIVO di -10.852001, indicando una sotto-espressione nei campioni di Heard 
# rispetto a quelli di Kidney. Il valore di p (PValue) è molto basso (4.350362e-09), 
# con un FDR altrettanto basso (0.0001004847), 
# suggerendo che questa differenza di espressione è altamente significativa statisticamente.

# Finally, a quick way to have an idea on how many DEG genes there are, 
# according to different FDR/FC thresholds:
summary(decideTests(qlfKB, p.value=0.05, adjust.method = "BH", lfc=0))
    # restituirà il numero di geni con un valore di p corretto < 0.05 senza il fold change (lfc = 0),
        
      #-1*Brain 1*kidney
#Down                4031. down-regulated in muscolo rispetto a cervello
#NotSig             15630. non mostrano una differenza significativa di espressione
#Up                  3377. up-regulated in muscolo rispetto a cervello


summary(decideTests(qlfKB, p.value=0.05, adjust.method = "BH", lfc=1))
    # restituirà il numero di geni con un valore di p corretto < 0.01 e un fold change ≥ 1.
summary(decideTests(qlfKB, p.value=0.01, adjust.method = "BH", lfc=0))
summary(decideTests(qlfKB, p.value=0.01, adjust.method = "BH",lfc=1))

summary(decideTests(qlfHB, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfHB, p.value=0.05, adjust.method = "BH", lfc=1))
summary(decideTests(qlfHB, p.value=0.01, adjust.method = "BH", lfc=0))
summary(decideTests(qlfHB, p.value=0.01, adjust.method = "BH",lfc=1))

summary(decideTests(qlfHK, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfHK, p.value=0.05, adjust.method = "BH", lfc=1))
summary(decideTests(qlfHK, p.value=0.01, adjust.method = "BH", lfc=0))
summary(decideTests(qlfHK, p.value=0.01, adjust.method = "BH",lfc=1))

### Now - a (small) challenge: build lists of genes that are up-regulated 
### in each tissue with respect the other two. That is - decide a threshold (e.g. FDR < 0.01) 
### and report all brain genes up-regulated in brain against kidney and heart. 

      # Filter genes that are up-regulated in brain vs kidney (negative logFC)
UP_brain_K<- topTags(qlfKB, n=Inf, adjust.method="BH", sort.by="PValue")
UP_brain_K <- UP_brain_K$table
UP_brain_K <- UP_brain_K[UP_brain_K$FDR < 0.01 & UP_brain_K$logFC < 0, ]
dim(UP_brain_K)
  
  
      # Filter genes that are up-regulated in brain vs heart (negative logFC)
UP_brain_H <- topTags(qlfHB, n=Inf, adjust.method="BH", sort.by="PValue")
UP_brain_H <- UP_brain_H$table
UP_brain_H <- UP_brain_H[UP_brain_H$FDR < 0.01 & UP_brain_H$logFC < 0, ]
dim(UP_brain_H)
      # Identify common genes up-regulated in brain in both comparisons
UP_brain <- intersect(rownames(UP_brain_K), rownames(UP_brain_H))
      # Report the list of genes
UP_brain_genes <- data.frame(Gene = UP_brain)


### Likewise for the other two tissues
UP_heart_K<- topTags(qlfHK, n=Inf, adjust.method="BH", sort.by="PValue")
UP_heart_K <- UP_heart_K$table
UP_heart_K <- UP_heart_K[UP_heart_K$FDR < 0.01 & UP_heart_K$logFC > 0, ]
dim(UP_heart_K)
UP_heart_B<- topTags(qlfHB, n=Inf, adjust.method="BH", sort.by="PValue")
UP_heart_B <- UP_heart_B$table
UP_heart_B <- UP_heart_B[UP_heart_B$FDR < 0.01 & UP_heart_B$logFC > 0, ]
dim(UP_heart_B)
UP_heart <- intersect(rownames(UP_heart_K), rownames(UP_heart_B))
UP_heart_genes <- data.frame(Gene = UP_heart)


UP_kidney_H<- topTags(qlfHK, n=Inf, adjust.method="BH", sort.by="PValue")
UP_kidney_H <- UP_kidney_H$table
UP_kidney_H <- UP_kidney_H[UP_kidney_H$FDR < 0.01 & UP_kidney_H$logFC < 0, ]
dim(UP_kidney_H)
UP_kidney_B<- topTags(qlfKB, n=Inf, adjust.method="BH", sort.by="PValue")
UP_kidney_B <- UP_kidney_B$table
UP_kidney_B <- UP_kidney_B[UP_kidney_B$FDR < 0.01 & UP_kidney_B$logFC > 0, ]
dim(UP_kidney_B)
UP_kidney <- intersect(rownames(UP_kidney_H), rownames(UP_kidney_B))
UP_kidney_genes <- data.frame(Gene = UP_kidney)

################################################################################
UP_brain_genes[1,]    ## [1] "CNTNAP4"
# For example, gene “CNTN2” is on the top of the DE lists, over-expressed in brain with respect to the other two tissues. 
# Let us plot the distribution of expression across the three complete datasets, as TPM:
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_heart)$TPM <- recount::getTPM(rse_heart)
assays(rse_kidney)$TPM <- recount::getTPM(rse_kidney)
which(rowData(rse_brain)$gene_name == "CNTNAP4")     #[1] 23734
#way we got the row number corresponding to gene “CNTNAP4”.
boxplot(assays(rse_brain)$TPM[23734,],
        assays(rse_heart)$TPM[23734,], 
        assays(rse_kidney)$TPM[23734,], outline=F, names=c("Brain", "heart", "kidney"), main = "CNTNAP4")

UP_heart_genes[1,]   ##[1] "TBX20"
which(rowData(rse_brain)$gene_name == "TBX20")    #[1] 45978
#we got the row number corresponding to gene “TBX20”.
boxplot(assays(rse_brain)$TPM[45978],
        assays(rse_heart)$TPM[45978,], 
        assays(rse_kidney)$TPM[45978,], outline=F, names=c("Brain", "heart", "kidney"), main = "TBX20")

UP_kidney_genes[1,]   ##[1] "KCNJ16"
which(rowData(rse_brain)$gene_name == "KCNJ16")    #[1] 25923
#we got the row number corresponding to gene “KCNJ16”.
boxplot(assays(rse_brain)$TPM[25923],
        assays(rse_heart)$TPM[25923,], 
        assays(rse_kidney)$TPM[25923,], outline=F, names=c("Brain", "heart", "kidney"), main = "KCNJ16")

# By default it outputs just the top 10 genes (the n parameter). 
# A very “quick and dirty” way to have the full table could be:
resultsHB <- topTags(qlfHB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsKB <- topTags(qlfKB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsKH <- topTags(qlfKH, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)




################################################################################

which(rowData(rse_brain)$gene_name == "FXYD2")    #[1] 14866
boxplot(assays(rse_brain)$TPM[14866],
        assays(rse_heart)$TPM[14866,], 
        assays(rse_kidney)$TPM[14866,], outline=F, names=c("Brain", "Heart", "Kidney"), main = "FXYD2")

# Finally, let us check whether one of the genes that were DE among these samples 
# are still DE if we consider the complete tissue datasets. 

################################################################################

# That can be saved in a file, e.g. for importing it in Excel:
write.table(resultsHB, "resultsHB.txt")
write.table(resultsKB, "resultsKB.txt")
write.table(resultsKH, "resultsKH.txt")

#link coi comandi del prof 
## http://159.149.160.56/GT_2022_GTEX/recount_DEG.html
## http://159.149.160.56/GT_2022_GTEX/recount3.html

UP_kidney_genes
UP_heart_genes
UP_brain_genes


brain_upregulated_first_500<-head(UP_brain_genes, 500)
heart_upregulated_first_500<-head(UP_heart_genes, 500)
kidney_upregulated_first_500<- head(UP_kidney_genes, 500)

write.csv(brain_upregulated_first_500, 
          file='./brain_upregulated_first_500.csv', row.names = FALSE)
file.exists('./brain_upregulated_first_500.csv')

write.csv(heart_upregulated_first_500, 
          file='./heart_upregulated_first_500.csv', row.names = FALSE)
write.csv(kidney_upregulated_first_500, 
          file='./kidney_upregulated_first_500.csv', row.names = FALSE)





