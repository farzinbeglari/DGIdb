# DGIdb
The drug gene interaction (PIK3CA)
library(httr)
library(jsonlite)
library(dplyr)
library(ggraph)
library(igraph)
library(tidyverse)

# Query DGIdb for PIK3CA
url <- 'https://dgidb.org/api/graphql'
query <- '
  query($genes: [String!]) {
    genes(names: $genes) {
      nodes {
        name
        interactions {
          drug {name approved}
          interactionScore
          interactionTypes {type}
        }
      }
    }
  }'

# Send the query to DGIdb
res <- POST(url,
            body = list(query = query, variables = list(genes = list('PIK3CA'))),
            encode = 'json')
data <- fromJSON(content(res, 'text', encoding = 'UTF-8'),
                 simplifyVector = F)

# Extract the results
rows <- list()
for (gene in data$data$genes$nodes) {
  for (i in gene$interactions) {
    rows[[length(rows) + 1]] <- data.frame(
      gene = gene$name,
      drug = i$drug$name,
      approved = isTRUE(i$drug$approved),
      score = ifelse(is.null(i$interactionScore), NA, i$interactionScore),
      types = paste(sapply(i$interactionTypes, `[[`, "type"), collapse = ", "),
      sources = paste(sapply(i$sources, `[[`, "sourceDbName"), collapse = ", "),
      pmids = paste(sapply(i$publications, `[[`, "pmid"), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }
}

df <- bind_rows(rows) %>%
  arrange(desc(score))

print(head(df, 20))
write.csv(df, "PIK3CA_DGIdb_interactions.csv", row.names = FALSE)

# Approved drugs visualization
df_app <- df %>%
  filter(approved) %>%
  separate_rows(types, sep = ',') %>%
  mutate(types = tolower(trimws(types))) %>%
  filter(!is.na(types), types != '', types != 'unknown')

table(df_app$types)

edges <- df_app %>%
  distinct(gene, drug, types) %>%
  rename(from = gene, to = drug)

nodes <- data.frame(name = unique(c(edges$from, edges$to))) %>%
  mutate(node_type = ifelse(name %in% edges$from, "Gene", "Drug"))

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

set.seed(42)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = types), width = 0.9, alpha = 0.85) +
  geom_node_point(aes(shape = node_type, fill = node_type),
                  size = 4, color = "black", stroke = 0.3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(Gene = 21, Drug = 22)) +
  scale_fill_manual(values = c(Gene = "#2E6F9E", Drug = "#C0392B")) +
  labs(title = "PIK3CA: approved drugs by known interaction type",
       edge_color = "Interaction type") +
  theme_void() +
  theme(legend.title = element_text(face = "bold"))

print(p)


PIK3CA is a key component of the insulin-mediated signaling pathway. Reduced expression may impair insulin signaling, potentially contributing to insulin resistance, abnormal glucose metabolism, and type 2 diabetes mellitus [1].

[1] Dong J, Zhang Y, Li YL, Wu LJ, Wang S, Chen N, Yan YX. Insulin signaling pathway related m6A methylated biomarker for type 2 diabetes and the potential modulation mechanism. Nutrition & Metabolism. 2026 Jan 31;23(1):28. doi: https://doi.org/10.1186/s12986-026-01086-4 
