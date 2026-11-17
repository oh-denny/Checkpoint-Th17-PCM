# Shared helper: Chord plot helpers; caller supplies plot_df and mac_col, neu_col, tar_col.
# Inputs: caller arguments/session objects; no automatic data loading.
# Outputs: functions and documented constants; plots only when called.
# Dependencies/prerequisites: see docs/SCRIPT_INVENTORY.md and SOFTWARE_REQUIREMENTS.md.

prepare_circos_data <- function(df, comparison_name, top_ligands_per_sender = 10, top_targets = 20, min_weight = 0.03) {
        x <- df %>% filter(comparison == comparison_name)

        lig_keep <- x %>%
        group_by(sender, ligand) %>%
        summarise(score = max(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
        group_by(sender) %>%
        slice_max(order_by = score, n = top_ligands_per_sender, with_ties = FALSE) %>%
        ungroup()

        x <- x %>% semi_join(lig_keep, by = c("sender", "ligand"))

        tgt_keep <- x %>%
        group_by(target) %>%
        summarise(score = max(weight, na.rm = TRUE), .groups = "drop") %>%
        slice_max(order_by = score, n = top_targets, with_ties = FALSE)

        x <- x %>% semi_join(tgt_keep, by = "target")

        x <- x %>% filter(weight >= min_weight)

        x <- x %>%
        mutate(from = ifelse(sender == "Macrophages", paste0("MAC__", ligand), paste0("NEU__", ligand)),
        to = paste0("T__", target),
        link_col = ifelse(sender == "Macrophages", adjustcolor(mac_col, alpha.f = 0.75), adjustcolor(neu_col, alpha.f = 0.75)))

        target_order <- x %>%
        group_by(target, to) %>%
        summarise(score = max(weight, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(score)) %>%
        pull(to)

        mac_order <- x %>%
        filter(sender == "Macrophages") %>%
        group_by(ligand, from) %>%
        summarise(score = max(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(score)) %>%
        pull(from)

        neu_order <- x %>%
        filter(sender == "Neutrophils") %>%
        group_by(ligand, from) %>%
        summarise(score = max(aupr_corrected, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(score)) %>%
        pull(from)

        sector_order <- c(target_order, mac_order, neu_order)

        grid.col <- c(setNames(rep(tar_col, length(target_order)), target_order),
        setNames(rep(mac_col, length(mac_order)), mac_order),
        setNames(rep(neu_col, length(neu_order)), neu_order))

        links <- x %>% select(from, to, value = weight, link_col)

        list(links = links, order = sector_order, grid.col = grid.col)
}



draw_circos_plot <- function(df, comparison_name, main_text, top_ligands_per_sender = 10, top_targets = 20, min_weight = 0.03) {
        circos.clear()

        cc <- prepare_circos_data(df, comparison_name, top_ligands_per_sender = top_ligands_per_sender, top_targets = top_targets, min_weight = min_weight)

        if(nrow(cc$links) == 0) stop(paste("No links found for", comparison_name))

        gap_vec <- c(rep(2, length(cc$order) - 1), 8)

        circos.par(start.degree = 90, gap.after = gap_vec, track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0))

        chordDiagram(as.data.frame(cc$links[, c("from", "to", "value")]),
            order = cc$order,
            grid.col = cc$grid.col,
            col = cc$links$link_col,
            transparency = 0.15,
            annotationTrack = "grid",
            preAllocateTracks = list(track.height = 0.14),
            link.sort = TRUE,
            link.largest.ontop = TRUE,
            directional = 0)

        circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
            sn <- get.cell.meta.data("sector.index")
            xl <- get.cell.meta.data("xlim")
            yl <- get.cell.meta.data("ylim")
            lab <- gsub("^MAC__", "", sn)
            lab <- gsub("^NEU__", "", lab)
            lab <- gsub("^T__", "", lab)
            circos.text(x = mean(xl), y = yl[1] + 0.1, labels = lab, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2, font = 3)
            })

        title(main = main_text, cex.main = 1.4, font.main = 2)
}
