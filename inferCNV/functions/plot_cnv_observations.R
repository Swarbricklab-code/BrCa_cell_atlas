infercnv:::.plot_cnv_observations
function (infercnv_obj, obs_data, col_pal, contig_colors, contig_labels,
    contig_names, contig_seps, num_obs_groups, file_base_name,
    output_filename_prefix, cnv_title, cnv_obs_title, contig_lab_size = 1,
    contigs, cluster_contig = NULL, obs_annotations_groups, obs_annotations_names,
    grouping_key_coln, cluster_by_groups, breaksList, x.center,
    hclust_method = "ward.D", testing = FALSE, layout_lmat = NULL,
    layout_lhei = NULL, layout_lwid = NULL)
{
    flog.info("plot_cnv_observation:Start")
    flog.info(paste("Observation data size: Cells=", nrow(obs_data),
        "Genes=", ncol(obs_data), sep = " "))
    observation_file_base <- paste(file_base_name, sprintf("%s.observations.txt",
        output_filename_prefix), sep = .Platform$file.sep)
    hcl_desc <- "General"
    hcl_group_indices <- 1:ncol(obs_data)
    if (!is.null(cluster_contig)) {
        hcl_contig_indices <- which(contig_names == cluster_contig)
        if (length(hcl_group_indices) > 0) {
            hcl_group_indices <- hcl_contig_indices
            hcl_desc <- cluster_contig
            flog.info(paste("plot_cnv_observation:Clustering only by contig ",
                cluster_contig))
        }
        else {
            flog.warn(paste("plot_cnv_observations: Not able to cluster by",
                cluster_contig, "Clustering by all genomic locations.",
                "To cluster by local genomic location next time",
                "select from:", unique(contig_names), collapse = ",",
                sep = " "))
        }
    }
    flog.info(paste("clustering observations via method: ", hclust_method,
        sep = ""))
    obs_dendrogram <- list()
    ordered_names <- NULL
    isfirst <- TRUE
    hcl_obs_annotations_groups <- vector()
    obs_seps <- c()
    sub_obs_seps <- c()
    if (!is.null(infercnv_obj@tumor_subclusters)) {
        for (i in seq(1, max(obs_annotations_groups))) {
            obs_dendrogram[[i]] = as.dendrogram(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])
            ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups ==
                i), hcl_group_indices])[(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])$order])
            obs_seps <- c(obs_seps, length(ordered_names))
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups,
                rep(i, length(which(obs_annotations_groups ==
                  i))))
            if (isfirst) {
                write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]),
                  file = paste(file_base_name, sprintf("%s.observations_dendrogram.txt",
                    output_filename_prefix), sep = .Platform$file.sep))
                isfirst <- FALSE
            }
            else {
                write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]],
                  file = paste(file_base_name, sprintf("%s.observations_dendrogram.txt",
                    output_filename_prefix), sep = .Platform$file.sep),
                  append = TRUE))
            }
        }
        if (length(obs_dendrogram) > 1) {
            obs_dendrogram <- do.call(merge, obs_dendrogram)
        }
        else {
            obs_dendrogram <- obs_dendrogram[[1]]
        }
        split_groups <- rep(1, dim(obs_data)[1])
        names(split_groups) <- ordered_names
        for (subtumor in infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]) {
            length(subtumor)
            sub_obs_seps <- c(sub_obs_seps, (sub_obs_seps[length(sub_obs_seps)] +
                length(subtumor)))
        }
    }
    else if (cluster_by_groups) {
        for (i in seq(1, max(obs_annotations_groups))) {
            gene_indices_in_group <- which(obs_annotations_groups ==
                i)
            num_genes_in_group <- length(gene_indices_in_group)
            flog.info(sprintf("Number of genes in group(%d) is %d",
                i, num_genes_in_group))
            if (num_genes_in_group < 2) {
                flog.info(sprintf("Skipping group: %d, since less than 2 entries",
                  i))
                next
            }
            data_to_cluster <- obs_data[gene_indices_in_group,
                hcl_group_indices, drop = FALSE]
            flog.info(paste("group size being clustered: ", paste(dim(data_to_cluster),
                collapse = ","), sep = " "))
            group_obs_hcl <- hclust(dist(data_to_cluster), method = hclust_method)
            ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups ==
                i), hcl_group_indices])[group_obs_hcl$order])
            if (isfirst) {
                write.tree(as.phylo(group_obs_hcl), file = paste(file_base_name,
                  sprintf("%s.observations_dendrogram.txt", output_filename_prefix),
                  sep = .Platform$file.sep))
                isfirst <- FALSE
            }
            else {
                write.tree(as.phylo(group_obs_hcl), file = paste(file_base_name,
                  sprintf("%s.observations_dendrogram.txt", output_filename_prefix),
                  sep = .Platform$file.sep), append = TRUE)
            }
            group_obs_dend <- as.dendrogram(group_obs_hcl)
            obs_dendrogram[[length(obs_dendrogram) + 1]] <- group_obs_dend
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups,
                rep(i, length(which(obs_annotations_groups ==
                  i))))
            obs_seps <- c(obs_seps, length(ordered_names))
        }
        if (length(obs_dendrogram) > 1) {
            obs_dendrogram <- do.call(merge, obs_dendrogram)
        }
        else {
            obs_dendrogram <- obs_dendrogram[[1]]
        }
        split_groups <- rep(1, dim(obs_data)[1])
        names(split_groups) <- ordered_names
        sub_obs_seps = obs_seps
    }
    else {
        obs_hcl <- hclust(dist(obs_data[, hcl_group_indices]),
            method = hclust_method)
        write.tree(as.phylo(obs_hcl), file = paste(file_base_name,
            sprintf("%s.observations_dendrogram.txt", output_filename_prefix),
            sep = .Platform$file.sep))
        obs_dendrogram <- as.dendrogram(obs_hcl)
        ordered_names <- row.names(obs_data)[obs_hcl$order]
        split_groups <- cutree(obs_hcl, k = num_obs_groups)
        split_groups <- split_groups[ordered_names]
        hcl_obs_annotations_groups <- obs_annotations_groups[obs_hcl$order]
        flog.info("plot_cnv_observation:Writing observations by grouping.")
        for (cut_group in unique(split_groups)) {
            group_memb <- names(split_groups)[which(split_groups ==
                cut_group)]
            memb_file <- file(paste(file_base_name, paste(hcl_desc,
                "HCL", cut_group, "members.txt", sep = "_"),
                sep = .Platform$file.sep))
            write.table(obs_data[group_memb, ], memb_file)
            ordered_memb <- which(ordered_names %in% group_memb)
            if (is.null(obs_seps)) {
                obs_seps <- c(length(ordered_memb))
            }
            else {
                obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] +
                  length(ordered_memb)))
            }
        }
        obs_seps <- c(obs_seps, length(ordered_names))
        sub_obs_seps = obs_seps
    }
    if (length(obs_seps) > 1) {
        obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) -
            1):1]
    }
    row_groupings <- get_group_color_palette()(length(table(split_groups)))[split_groups]
    row_groupings <- cbind(row_groupings, get_group_color_palette()(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])
    annotations_legend <- cbind(obs_annotations_names, get_group_color_palette()(length(table(hcl_obs_annotations_groups))))
    flog.info("plot_cnv_observation:Writing observation groupings/color.")
    groups_file_name <- file.path(file_base_name, sprintf("%s.observation_groupings.txt",
        output_filename_prefix))
    file_groups <- cbind(split_groups, row_groupings[, 1], hcl_obs_annotations_groups,
        row_groupings[, 2])
    colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color",
        "Annotation Group", "Annotation Color")
    write.table(file_groups, groups_file_name)
    contigSepList <- create_sep_list(row_count = nrow(obs_data),
        col_count = ncol(obs_data), row_seps = obs_seps, col_seps = contig_seps)
    obs_data <- obs_data[ordered_names, ]
    orig_row_names <- row.names(obs_data)
    row.names(obs_data) <- rep("", nrow(obs_data))
    heatmap_thresholds_file_name <- file.path(file_base_name,
        sprintf("%s.heatmap_thresholds.txt", output_filename_prefix))
    write.table(breaksList, heatmap_thresholds_file_name, row.names = FALSE,
        col.names = FALSE)
    data_observations <- heatmap.cnv(obs_data, Rowv = obs_dendrogram,
        Colv = FALSE, cluster.by.row = TRUE, cluster.by.col = FALSE,
        main = cnv_title, ylab = cnv_obs_title, margin.for.labCol = 2,
        xlab = "Genomic Region", key = TRUE, labCol = contig_labels,
        cexCol = contig_lab_size, notecol = "black", density.info = "histogram",
        denscol = "blue", trace = "none", dendrogram = "row",
        cexRow = 0.8, breaks = breaksList, scale = "none", x.center = x.center,
        color.FUN = col_pal, if.plot = !testing, sepList = contigSepList,
        sep.color = c("black", "black"), sep.lty = 1, sep.lwd = 1,
        RowIndividualColors = row_groupings, annotations_legend = annotations_legend,
        grouping_key_coln = grouping_key_coln, ColIndividualColors = contig_colors,
        key.title = "Distribution of Expression", key.xlab = "Modified Expression",
        key.ylab = "Count", force_lmat = layout_lmat, force_lwid = layout_lwid,
        force_lhei = layout_lhei)
    if (class(obs_data) %in% c("matrix", "data.frame")) {
        flog.info(paste("plot_cnv_references:Writing observation data to",
            observation_file_base, sep = " "))
        row.names(obs_data) <- orig_row_names
        write.table(t(obs_data[data_observations$rowInd, data_observations$colInd]),
            file = observation_file_base)
    }
}
<bytecode: 0x1aa9af18>
<environment: namespace:infercnv>