### [ WHATSTHATCELL FUNCTIONS ] #####
select_initial_cells <- function(df_expression,
                                 marker_dict,
                                 number_cells = 2) {
  cell_types <- names(marker_dict)
  unannotated_cells <- df_expression %>% filter(is.na(cell_type))
  annotated_cells <- df_expression %>% filter(!is.na(cell_type))
  left_cell_types <-
    setdiff(cell_types, annotated_cells$cell_type)
  
  selected_cells <- lapply(left_cell_types, function(x) {
    #Negative markers
    neg_markers <- marker_dict[[x]]$negative
    if (!is.null(neg_markers)) {
      neg_marker_expr <- unannotated_cells[, neg_markers]
      neg_avg_expression <-
        apply(neg_marker_expr, 1, calculate_expression_average)
    } else{
      neg_avg_expression <- 0
    }
    #Positive markers
    pos_markers <- marker_dict[[x]]$positive
    if (!is.null(pos_markers)) {
      pos_marker_expr <- unannotated_cells[, pos_markers]
      pos_avg_expression <-
        apply(pos_marker_expr, 1, calculate_expression_average)
    } else{
      pos_avg_expression <- 0
    }
    
    average_expression <- pos_avg_expression - neg_avg_expression
    unannotated_cells$average_expression <- average_expression
    unannotated_cells[order(unannotated_cells$average_expression, decreasing = T), ][1:number_cells, ]$X1
  }) %>% unlist()
  return(unique(selected_cells))
}

calculate_expression_average <- function(expressions) {
  return(sum(expressions) / length(expressions))
}

load_scs <- function(scs_expression) {
  df_expression <- as.data.frame(as.matrix(t(assays(scs_expression)$logcounts)))
  setDT(df_expression, keep.rownames = "X1")[]
  df_expression <- as_tibble(df_expression)
  return(df_expression)
}

calculate_entropy <- function(predictions) {
  entropy <- -(lapply(predictions, function(x) {
    if (x == 0) {
      return(0)
    } else {
      return(x * log(x, 2))
    }
  }) %>% unlist() %>% sum())
  return(entropy)
}

select_cells_classifier <- function(df_expression, markers, amount = 10) {
  annotated_cells <- df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    filter(cell_type != "Skipped", cell_type != "Unclear")
  
  left_cells <- df_expression %>% 
    filter(is.na(cell_type))
  
  multiNomModelFit <- train(cell_type ~ ., 
                            data = annotated_cells[, c(markers, "cell_type")], 
                            method = "multinom")
  
  predicted_scores <- predict(multiNomModelFit, 
                              left_cells[, markers], type = "prob")
  
  entropies <- apply(predicted_scores, 1, calculate_entropy)
  
  left_cells$entropy <- entropies
  
  entropy_table <- data.frame(cbind(cell_id = left_cells$X1, 
                                    entropy = left_cells$entropy, 
                                    no_cells_annotated = nrow(annotated_cells)))
  
  selected_cells <- c(left_cells[order(left_cells$entropy, 
                                       decreasing = T),][1:amount, ]$X1)
  
  return_list <- list(selected_cells = selected_cells, 
                      entropy_table = entropy_table)
  return(return_list)
}


acc_wrap <- function(tt) {
  cell_types <- unique(union(tt$predicted_cell_type, tt$annotated_cell_type))
  
  tt$annotated_cell_type <- factor(tt$annotated_cell_type, levels = cell_types)
  tt$predicted_cell_type <- factor(tt$predicted_cell_type, levels = cell_types)
  
  bind_rows(
    tryCatch({kap(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({mcc(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({sensitivity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({specificity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL)
  )
}


cell_type_colours <- function() {
  pal <- c("#8B5B42", "#AF4EA9", "#FFB60A", "#0AC694", "#0024DD", "#6CC1FF",
           "#0496FF", "#1DA05B", "#E11E00", "#A78882", "#BD93D8", "#fff53d")
  
  celltype_colours <- c(
    "B cell" = pal[2],
    "Cytotoxic T cell" = pal[3],
    "CD4+ T cell" = pal[12],
    "CD16+ monocyte" = pal[4],
    "Dendritic cell" = pal[5],
    "Plasmacytoid dendritic cell" = pal[6],
    "CD14+ monocyte" = pal[8], 
    "Megakaryocyte" = pal[9],
    "Natural killer cell" = pal[10],
    'unassigned' = "grey60"
  )
  celltype_colours
}


whatsthatcell_theme <- function(){
  theme_bw()
}
