app <- ShinyDriver$new("../../")
app$snapshotInit("query_genes")

app$setInputs(menu = "Query genes")
app$setInputs(select_gene = "SERPINA3")
app$setInputs(select_gene = c("A1BG", "SERPINA3"))
app$setInputs(select_gene = c("A1BG", "A1CF", "SERPINA3"))
app$snapshot()
