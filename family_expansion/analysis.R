# Dependencies
library(Biostrings)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)

# Load functions 
#`blast_homologs_nt`available from: https://github.com/jdieramon/my_scripts/blob/master/genomics/genomics.R#L70

blast_homologs_nt <- function(hitFile, ident, cover){
  
  ### ''' Take a csv table with Blast hits "Align two or more sequences"
  ### return subset filtered by %Identity & %Coverage (longer sequence)
  # Ex. Duplication Analysis
  
  
  # read downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # change colnames
  # en blast-2 sequences contra ellas mismas no hay 13 columnas sino 12 (-"%positives")
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score")
  
  # make a table_length for queries & subjects 
  tlengths <- hit %>% 
    filter(query == subject) %>% 
    select(query, align_length) %>% 
    arrange(query, desc(align_length)) %>% 
    distinct(query, .keep_all = TRUE) %>% 
    dplyr::rename(length = align_length)
  
  
  # combine hit_table with queries & subjects lengths 
  # estimate coverage (longer sequence)
  as_tibble(hit) %>% 
    select(-c(align_length, bit_score, mismatches, gap_opens)) %>% 
    filter(query != subject) %>% 
    right_join(tlengths, by = "query") %>% 
    dplyr::rename(len.query = length) %>% 
    right_join(tlengths %>% dplyr::rename(subject = query), by = "subject") %>% 
    dplyr::rename(len.sub = length) %>% 
    mutate(cov.query = (q.end - q.start + 1)/len.query, 
           cov.sub = (s.end - s.start + 1)/len.sub) %>% 
    mutate(cov_long = case_when(len.query > len.sub ~ cov.query*100, 
                                len.query < len.sub ~ cov.sub*100, 
                                TRUE ~ cov.query*100)) %>% 
    select(-c(q.start, q.end, s.start, s.end, len.query, len.sub, cov.query, cov.sub)) %>% 
    filter(identity >= ident, cov_long >= cover) %>% 
    select(query, subject, identity, cov_long, evalue) %>% 
    arrange(identity)
  
}


# Load data (fasta file)
cds <- readDNAStringSet("myCDS.fasta")

# Tidy names 
names(cds) <- str_sub(names(cds), start = 1, end = 12)

# re-write file with tidy names 
cds 
writeXStringSet(cds, filepath="CDStidynames.fasta", format="fasta") 


# run BLAST pair-wise 
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# download Hit table (.csv)
hitFile <- "3BJET9JM114-Alignment-HitTable.csv"

df <- blast_homologs_nt(hitFile, ident = 70, cover = 70)




# problema clásico de detección de componentes conexos en un grafo. En este caso:
#   
#   Cada secuencia (query o subject) es un nodo.
# 
#   Cada fila representa una conexión (arista) entre dos nodos.
# 
#   Una familia corresponde a un componente conexo del grafo.



# Crear un grafo no dirigido con los pares query-subject
g <- graph_from_data_frame(df %>% select(query, subject), directed = FALSE)

# Detectar componentes conexos
componentes <- components(g)

# Extraer las familias
familias <- split(names(componentes$membership), componentes$membership)

# Asignar nombre a las familias
names(familias) <- paste0("familia_", seq_along(familias))

# Resultado: lista con las familias y sus miembros
print(familias)



# Visualización básica
plot(
  g,
  vertex.color = componentes$membership + 1,  # Colores por familia
  vertex.label = V(g)$name,                  # Etiquetas con los nombres
  vertex.size = 30,
  edge.arrow.size = 0.5,
  layout = layout_with_fr,                   # Algoritmo de layout
  main = "Familias de secuencias (componentes conexos)"
)


# Visualización ggraph 
# Convertir a objeto tidygraph
tg <- as_tbl_graph(g)

# Añadir la familia como atributo
tg <- tg %>% mutate(familia = factor(componentes$membership))

# Graficar con ggraph
ggraph(tg, layout = "fr") +
  geom_edge_link(linetype = "dashed", edge_colour = "gray40", alpha = 0.6) +
  geom_node_point(aes(color = familia), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_minimal() + xlab("") + ylab("") +
  labs(title = "Familias de secuencias (componentes conexos)")

# layout = "fr" usa el algoritmo Fruchterman-Reingold, que distribuye los nodos de forma estética y clara. Las posiciones de los nodos (X/Y) se calculan automáticamente mediante a
# lgoritmos de layout como Fruchterman-Reingold o Kamada-Kawai, que solo intentan distribuir los nodos de forma legible. Por tanto, los valores de los ejes son arbitrarios y no representan 
# variables medibles y no tienen significado geométrico real.   





# Análisis thresholds 80%-90% 
df <- blast_homologs_nt(hitFile, ident = 80, cover = 80)
df <- blast_homologs_nt(hitFile, ident = 90, cover = 90)





