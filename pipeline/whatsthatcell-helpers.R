### [ WHATSTHATCELL FUNCTIONS ] #####

### [ ACTIVE LEARNING ] #####
select_initial_cells <- function(df_expression,
                                 marker_dict,
                                 number_cells = 20) {
  
  cell_types <- names(marker_dict)
  
  selected_cells <- lapply(cell_types, function(x) {
    #Negative markers
    neg_markers <- marker_dict[[x]]$negative
    if (!is.null(neg_markers)) {
      neg_marker_expr <- df_expression[, neg_markers]
      neg_avg_expression <-
        apply(neg_marker_expr, 1, calculate_expression_average)
    } else{
      neg_avg_expression <- 0
    }
    #Positive markers
    pos_markers <- marker_dict[[x]]$positive
    if (!is.null(pos_markers)) {
      pos_marker_expr <- df_expression[, pos_markers]
      pos_avg_expression <-
        apply(pos_marker_expr, 1, calculate_expression_average)
    } else{
      pos_avg_expression <- 0
    }
    
    average_expression <- pos_avg_expression - neg_avg_expression
    df <- tibble(marker = average_expression)
    colnames(df) <- names(marker_dict[x])
    df
  }) %>% bind_cols()
  
  selected_cells <- selected_cells |> 
    rowwise() |> 
    mutate(max_mean_expression = max(across()),
           max_cell_type = cell_types[which.max(c_across(everything(selected_cells)))]) |> 
    ungroup() |> 
    mutate(X1 = df_expression$X1) |> 
    select(X1, max_cell_type, max_mean_expression) |> 
    arrange(-max_mean_expression)
    
  cell_types_to_sample <- recycle(cell_types, number_cells)
  
  lapply(cell_types_to_sample, function(i){
    cell <- filter(selected_cells, max_cell_type == i) |> 
      slice_head(n = 1) |> 
      pull(X1)
    
    selected_cells <<- filter(selected_cells, !(X1 %in% cell))
    cell
  }) |> unlist()
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

# entropy based selection of cells
select_entropy <- function(entropy_df, method, amount, random_selection){
  # Determine how many cells should be selected non-randomly
  select_amount <- round(amount * (1-random_selection), 0)
  
  # Determine how many cells should be selected randomly
  random_amount <- amount - select_amount
  
  # Select the cells with the highest entropy
  if(method == "highest_entropy"){
    cells <- c(entropy_df[order(entropy_df$entropy, 
                                decreasing = T),][1:select_amount, ]$X1)
  }
  
  # Find appropriate quantile threshold
  else if(method == "0.95_quant_entropy"){
    quantile <- quantile(entropy_df$entropy, 0.95)
  }else if(method == "0.75_quant_entropy"){
    quantile <- quantile(entropy_df$entropy, 0.75)
  }
  
  # Filter cells based on quantiles
  if(grepl("_quant_entropy", method)){
    ordered_cells <- select(entropy_df, X1, entropy) %>% 
      arrange(entropy) %>% 
      filter(entropy > quantile)
    
    if(nrow(ordered_cells) > amount){
      cells <- ordered_cells$X1[1:select_amount]
    }else{
      cells <- ordered_cells$X1
    }
  }
  
  # Randomly select additional cells
  if(random_selection != 0){
    if(nrow(entropy_df) > random_amount){
      random_cells <- sample(entropy_df$X1, random_amount)
    }else{
      random_cells <- entropy_df$X1
    }
    cells <- c(cells, random_cells)
  }
  
  cells
}

# maximum probability based selection of cells
select_maxp <- function(maxp_df, method = "lowest_maxp", amount, random_selection){
  # Determine how many cells should be selected non-randomly
  select_amount <- round(amount * (1-random_selection), 0)
  
  # Determine how many cells should be selected randomly
  random_amount <- amount - select_amount
  
  # Select the cells with the highest maxp
  if(method == "lowest_maxp"){
    cells <- maxp_df[order(maxp_df$maxp, decreasing = F),][1:select_amount, ]$X1
  }
  
  
  # Find appropriate quantile threshold
  else if(method == "0.05_quant_maxp"){
    quantile <- quantile(maxp_df$maxp, 0.05)
  }else if(method == "0.25_quant_maxp"){
    quantile <- quantile(maxp_df$maxp, 0.25)
  }
  
  # Filter cells based on quantiles
  if(grepl("_quant_maxp", method)){
    ordered_cells <- select(maxp_df, X1, maxp) %>%
      arrange(maxp) %>%
      filter(maxp > quantile)
    
    if(nrow(ordered_cells) > amount){
      cells <- ordered_cells$X1[1:select_amount]
    }else{
      cells <- ordered_cells$X1
    }
  }
  
  # Randomly select additional cells
  if(random_selection != 0){
    if(nrow(maxp_df) > random_amount){
      random_cells <- sample(maxp_df$X1, random_amount)
    }else{
      random_cells <- maxp_df$X1
    }
    cells <- c(cells, random_cells)
  }
  
  cells
}

select_cells_classifier <- function(df_expression, AL_method, selection_method, amount = 10, 
                                    random_selection, selection_criterion = "entropy") {
  annotated_cells <- df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    filter(cell_type != "Skipped", cell_type != "Unclear")
  
  left_cells <- df_expression %>% 
    filter(is.na(cell_type))
  
  training_set_size <- nrow(annotated_cells)
  if(training_set_size < 50){
    bootstrap_reps <- 1000
  }else if(training_set_size >=50 | training_set_size < 100){
    bootstrap_reps <- 50
  }else{
    bootstrap_reps <- 25
  }
  tctrl <- trainControl(method = "boot", number = bootstrap_reps)
  ModelFit <- train(cell_type ~ ., 
                    data = select(annotated_cells, -X1),
                    method = AL_method,
                    trace = FALSE,
                    trControl = tctrl)
  
  predicted_scores <- predict(ModelFit, 
                              select(left_cells, -X1, -cell_type),
                              type = "prob")
  
  if(selection_criterion == "entropy"){
    entropies <- apply(predicted_scores, 1, calculate_entropy)
    
    left_cells$entropy <- entropies
    
    selected_cells <- select_entropy(entropy_df = left_cells, 
                                     method = selection_method, 
                                     amount = amount, 
                                     random_selection = random_selection)
    
    criterion_table <- data.frame(cell_id = left_cells$X1, 
                                  criterion_val = left_cells$entropy, 
                                  no_cells_annotated = nrow(annotated_cells),
                                  criterion = selection_criterion)
    
  }else if(selection_criterion == "maxp"){
    max_p_idx <- apply(predicted_scores, 1, which.max)
    max_p <- lapply(seq_len(nrow(predicted_scores)), function(x){
      predicted_scores[x, max_p_idx[x]]
    }) %>% unlist()
    
    left_cells$maxp <- max_p
    
    selected_cells <- select_maxp(maxp_df = left_cells, 
                                  method = selection_method, 
                                  amount = amount, 
                                  random_selection = random_selection)
    
    criterion_table <- data.frame(cell_id = left_cells$X1, 
                                  criterion_val = left_cells$maxp, 
                                  no_cells_annotated = nrow(annotated_cells),
                                  criterion = selection_criterion)
  }
  
  
  return_list <- list(selected_cells = selected_cells, 
                      criterion_table = criterion_table)
  return(return_list)
}


cell_ranking_wrapper <- function(df, markers, number_cells = 20){
  # Get initial set of cells based on their marker expression ranking
  ranked_cells <- select_initial_cells(df, markers$cell_types, number_cells)
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(ranked_cells, df$X1)
  df$cell_type[to_assign_index] <- df$gt_cell_type[to_assign_index]
  df$iteration <- 0
  
  df
}


active_learning_wrapper <- function(df, AL_method, selection_method, iteration, 
                                    entropies = NULL, random_selection, selection_criterion = "entropy"){
  # AL selected cells
  AL <- select_cells_classifier(df_expression = df, 
                                AL_method = AL_method, 
                                selection_method = selection_method,
                                amount = 10, 
                                random_selection = random_selection,
                                selection_criterion = selection_criterion)
  
  list(new_cells = AL$selected_cells, criterion_table = AL$criterion_table)
}

corrupt_labels <- function(df_expression, all_cell_types, prop_corrupt){
  num_corrupt <- round((nrow(df_expression) * as.numeric(prop_corrupt)), 0)
  corrupt_idx <- sample(1:nrow(df_expression), num_corrupt)

  for(i in 1:length(all_cell_types)){
    tc <- all_cell_types[i]

    to_corrupt <- which(df_expression$gt_cell_type == tc)
    to_corrupt <- intersect(to_corrupt, corrupt_idx)
    if(length(to_corrupt) > 0){
      subset_cell_types <- all_cell_types[!(tc == all_cell_types)]
      df_expression$corrupted[to_corrupt] <- sample(subset_cell_types, 
                                                    length(to_corrupt), 
                                                    replace = TRUE)
    }
  }
  
  df_expression$gt_cell_type[!is.na(df_expression$corrupted)] <- 
    df_expression$corrupted[!is.na(df_expression$corrupted)]
  
  df_expression
}


### [ ACCURACIES ] #####
acc_wrap <- function(tt) {
  cell_types <- unique(union(tt$predicted_cell_type, tt$annotated_cell_type))
  
  tt$annotated_cell_type <- factor(tt$annotated_cell_type, levels = cell_types)
  tt$predicted_cell_type <- factor(tt$predicted_cell_type, levels = cell_types)
  
  bind_rows(
    tryCatch({kap(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({mcc(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({sensitivity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({f_meas(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
  )
}


### [ PLOTTING ] ####
createHeatmap <- function(sce,
                          cell_type_column = "cell_type",
                          assay = "logcounts",
                          thresh = 2,
                          title = NULL,
                          modality) {
  
  lc <- t(as.matrix(assay(sce, assay)))
  lc <- scale(lc)
  keep_cols <- apply(lc, 2, function(x) all(is.na(x) == FALSE)) %>% as.vector
  lc <- lc[, keep_cols]
  
  lc[lc > thresh] <- thresh
  lc[lc < -thresh] <- -thresh
  
  cell_types = colData(sce)[[ cell_type_column ]]
  
  celltype_annot <- HeatmapAnnotation(`Cell type` = cell_types, 
                                      which="column",
                                      col = list(`Cell type` = cell_type_colours(modality)))  
  
  if(is.null(title)){
    title <- title
  }else{
    title <- paste(title, "\nCell")
  }
  
  type_exprs <- Heatmap(t(lc), 
                        name = "Expression",
                        column_title = title,
                        col=viridis(100),
                        top_annotation = celltype_annot,
                        show_column_names = FALSE,
                        column_order = order(cell_types))
  type_exprs
}



cell_type_colours <- function(modality) {
  pal <- c("#8B5B42", "#AF4EA9", "#FFB60A", "#79E200", "#0024DD", "#6CC1FF",
           "#0496FF", "#1DA05B", "#E11E00", "#A78882", "#BD93D8", "#fff53d", '#FD4FBD')
  
  scRNASeq_colours <- c(
    "B cell" = pal[2],
    "Cytotoxic T cell" = pal[3],
    "CD4+ T cell" = pal[12],
    "CD16+ monocyte" = pal[4],
    "Dendritic cell" = pal[5],
    "CD14+ monocyte" = pal[8], 
    "Megakaryocyte" = pal[9],
    "Natural killer cell" = pal[10],
    "Plasmacytoid dendritic cell" = pal[6],
    'unassigned' = "grey60"
  )
  snRNASeq_colours <- c(
    "Fibroblast" = pal[1],
    "Tumor" = pal[2],
    "Endothelial" = pal[9],
    "Immune" = pal[3],
    "Ductal" = pal[5],
    "SmoothMuscle" = pal[10],
    "Schwann" = pal[4],
    "Atypical_Ductal" = pal[7],
    "Endocrine" = pal[8],
    "Acinar" = pal[6],
    'unassigned' = "grey60"
  )
  CyTOF_colours <- c(
    "B-cell Frac A-C (pro-B cells)" = pal[2],
    "IgD- IgMpos B cells" = pal[11],
    "IgM- IgD- B-cells" = pal[13],
    "CMP" = pal[3],
    'GMP' = pal[4],
    "CLP" = pal[8],
    "MPP" = pal[9],
    'unassigned' = "grey60"
  )
  
  if(modality == "scRNASeq"){
    scRNASeq_colours
  }else if(modality == "snRNASeq"){
    snRNASeq_colours
  }else if(modality == "CyTOF"){
    CyTOF_colours
  }
}


whatsthatcell_theme <- function(){
  theme_bw()
}



