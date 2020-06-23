# SERVER
server = function(input, output, session) {
  
#### Query genes ####
# boxplots
output$gene_regulation_boxplot = renderPlotly({
  if (!is.null(input$select_gene)) {
    gene_regulation_boxplot = contrasts %>%
      # filter(gene %in% c("MYC")) %>%
      filter(gene %in% input$select_gene) %>%
      ggplot(aes(x = fct_reorder(gene, t, median), y = t, label = study)) +
      geom_boxplot() +
      geom_point(aes(color = study), size = 4, alpha = 0.6) +
      theme_classic() +
      labs(x = "Queried genes", y = "t-value", color = "Study") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size= 11)) +
      geom_hline(yintercept = 0, color = "grey", linetype = 2) +
      # geom_vline(xintercept = length(HFgenes_up)+0.5, color = "black", linetype =1) +
      ggtitle("")
    
    ggplotly(gene_regulation_boxplot, tooltip = c("label"))

  }
})
  
# plot for rank positions
output$rank_position = renderPlotly({
  if (!is.null(input$select_gene)) {
    sub_ranks = ranks %>%
      filter(gene %in% input$select_gene)
    
    max_rank = max(ranks$rank)
    
    rank_plot = sub_ranks %>%
      ggplot(aes(x=rank, y=1, label = gene)) +
      geom_segment(mapping = aes(y = 0.5, xend = rank, yend = 1.5), 
                   size = 0.5, color=aachen_color("red"), alpha=0.5) +
      geom_hline(yintercept = 1, size=1.5) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(),
        axis.line.x = element_blank()
      ) +
      labs(x = "Rank", y=NULL) +
      xlim(1, max_rank)
    
    ggplotly(rank_plot, tooltip = c("label", "x"))
  }
})
 
# distribution of mean t-values 
output$mean_t_dist = renderPlotly({
  if (!is.null(input$select_gene)) {
    # density
    sub_ranks = ranks %>%
      filter(gene %in% input$select_gene)
    dens = ranks %>%
      ggplot(aes(x=mean_t, label = gene)) +
      stat_density(geom = "line") +
      geom_rug(data = sub_ranks, color=aachen_color("red"), size=1) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "mean t-value", y="density")
    
    ggplotly(dens, tooltip = c("label"))
  }
})
  
# subsetted data frame of gene expression
output$individual_sub = DT::renderDataTable({
  if (!is.null(input$select_gene)) {
    contrasts %>%
      filter(gene %in% input$select_gene) %>%
      mutate(AveExpr = signif(AveExpr,3),
             logFC = signif(logFC,3),
             t = signif(t,3),
             B = signif(B,3),
             P.Value = scientific(P.Value),
             adj.P.Val = scientific(adj.P.Val),
             study = as_factor(study)) %>%
      DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                    extensions = "Buttons", rownames = F,
                    option = list(scrollX = T, 
                                  autoWidth = T, 
                                  dom = "Bfrtip",
                                  buttons = c("copy", "csv", "excel")))
  }
})
  
# subsetted dataframe of meta analysis
output$summary_sub = DT::renderDataTable({
  if (!is.null(input$select_gene)) {
    ranks %>%
      filter(gene %in% input$select_gene) %>%
      mutate(mean_lfc = signif(mean_lfc,3),
             mean_t = signif(mean_t,3),
             fisher_pvalue = scientific(fisher_pvalue)) %>%
      DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                    extensions = "Buttons", rownames = F,
                    option = list(scrollX = T, 
                                  autoWidth = T, 
                                  dom = "Bfrtip",
                                  buttons = c("copy", "csv", "excel")))
  }
})

#### Input data ####
# load user input (gene sets)
gs = reactive({
  if (input$take_example_data == F) {
    shinyjs::enable("user_input")
    inFile = input$user_input
    if (is.null(inFile)){
      return(NULL)
    }
    read_csv(inFile$datapath)
  } else {
    shinyjs::disable("user_input")
    example_geneset
  }
})

# get which signature should be used (directed vs undirected)
signature = reactive({
  switch(input$signature_source,
         "directed" = directed_signature,
         "undirected" = undirected_signature)
})

# perform GSEA with plotting
gsea_res = eventReactive(input$submit, {
  if (ncol(gs()) == 1) {
    res = fgsea(list(geneset = gs()$gene), deframe(signature()), nperm = 1000)
    p = make_gsea_plot(signature(), gs(), weight)
  } else if (ncol(gs()) == 2) {
    list_gs = gs() %>%
      split(.$geneset) %>%
      map(pull, gene)
    res = fgsea(list_gs, deframe(signature()), nperm = 1000)
    
    plot_df = gs() %>%
      nest(set = gene) %>%
      mutate(plot = pmap(., .f = function(geneset, set, ...) {
        make_gsea_plot(signature(), set, weight)
      }))
    
    p = plot_grid(plotlist = plot_df$plot, labels = plot_df$geneset)
  }
  
  df = res %>% 
    rename(geneset = pathway) %>%
    as_tibble() %>%
    select(-leadingEdge) %>%
    mutate(signature = input$signature_source)
  
  list(df = df, p = p)
})

# gsea results as table
output$gsea_res_table = DT::renderDataTable({
  gsea_res()$df %>%
    mutate(NES = signif(NES, 3),
           ES = signif(ES, 3),
           pval = scientific(pval),
           padj = scientific(padj)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

# gsea results as plots
output$gsea_res_plots = renderPlot({
  gsea_res()$p
})


#### Meta analysis results ####
output$individual = DT::renderDataTable({
  contrasts %>%
    mutate(AveExpr = signif(AveExpr,3),
           logFC = signif(logFC,3),
           t = signif(t,3),
           B = signif(B,3),
           P.Value = scientific(P.Value),
           adj.P.Val = scientific(adj.P.Val),
           study = as_factor(study)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$summary = DT::renderDataTable({
  ranks %>%
    mutate(mean_lfc = signif(mean_lfc,3),
           mean_t = signif(mean_t,3),
           fisher_pvalue = scientific(fisher_pvalue)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
  })

#### Functional analysis ####
# progeny
output$progeny_table = DT::renderDataTable({
  progeny %>%
    mutate(progeny_scores = signif(progeny_scores, 3),
           progeny_pvals = scientific(progeny_pvals),
           pathway = factor(pathway, levels = sort(pathway))) %>%
    select(pathway, progeny_scores, progeny_pvals) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

# dorothea
output$dorothea_table = DT::renderDataTable({
  dorothea %>%
    mutate(NES = signif(NES, 3),
           pvalue = scientific(pvalue),
           adj_pvalue = scientific(adj_pvalue),
           RegulonName = factor(RegulonName, levels = sort(RegulonName))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

# gsea
output$gsea_table = DT::renderDataTable({
  gsea %>%
    mutate(NES = signif(NES, 3),
           ES = signif(ES, 3),
           pval = scientific(pval),
           padj = scientific(padj),
           database = factor(database, levels = sort(unique(database))),
           pathway = factor(pathway, levels = sort(unique(pathway)))) %>%
    select(-leadingEdge) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

# gsea with microRNA
output$mi_gsea_table = DT::renderDataTable({
  mi_gsea %>%
    mutate(NES = signif(NES, 3),
           pvalue = scientific(pvalue),
           adj_pvalue = scientific(adj_pvalue),
           RegulonName = factor(RegulonName, 
                                levels = sort(unique(RegulonName)))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})


#### Study overview ####
output$overview = DT::renderDataTable({
  overview %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel"),
                                pageLength = 14))
})

disable("submit")

observeEvent({
  input$user_input
  input$take_example_data
  }, {
    if (input$take_example_data == T) {
      enable("submit")
      disable("user_input")
    } else if (input$take_example_data == F & is.null(input$user_input)) {
      disable("submit")
      enable("user_input")
    } else if (input$take_example_data == F & !is.null(input$user_input))
      enable("submit")
  })

hide("loading-content", TRUE, "fade")  

}


