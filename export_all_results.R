# =============================================================================
# EXPORT ALL RESULTS
# 把所有分析结果（表格 + 图）一次性导出到工作目录
# 在你的分析脚本末尾运行这段代码
# =============================================================================

library(writexl)
library(ggplot2)

# ── 1. 确认工作目录 ────────────────────────────────────────────────────────────
message("Working directory: ", getwd())
message("All outputs will be saved here.")


# ── 2. 导出所有结果表格到一个 Excel 文件 ──────────────────────────────────────
# 自动检测当前环境中存在哪些结果数据框

export_list <- list()

# BBD truncating results
if (exists("res_bbd"))     export_list[["BBD_overall"]]          <- fmt_table(res_bbd)
if (exists("res_nonprol")) export_list[["BBD_nonproliferative"]] <- fmt_table(res_nonprol)
if (exists("res_prol"))    export_list[["BBD_proliferative"]]    <- fmt_table(res_prol)

# DCIS / LCIS truncating results (if run)
if (exists("res_dcis"))    export_list[["DCIS_truncating"]]      <- fmt_table(res_dcis)
if (exists("res_lcis"))    export_list[["LCIS_truncating"]]      <- fmt_table(res_lcis)

# Missense results (if run)
if (exists("res_miss_bbd"))  export_list[["BBD_missense"]]       <- fmt_table(res_miss_bbd)
if (exists("res_miss_dcis")) export_list[["DCIS_missense"]]      <- fmt_table(res_miss_dcis)
if (exists("res_miss_lcis")) export_list[["LCIS_missense"]]      <- fmt_table(res_miss_lcis)

# Objective 2 hormonal factors (if run)
if (exists("results_obj2_df")) export_list[["Obj2_hormonal"]]    <- results_obj2_df

# Combined raw results
if (exists("results_all")) export_list[["All_raw"]]              <- results_all

# Write Excel
xlsx_path <- "dissertation_all_results.xlsx"
write_xlsx(export_list, path = xlsx_path)
message(sprintf("\n✓ Excel exported: %s  (%d sheets)", xlsx_path, length(export_list)))


# ── 3. 重新生成并保存所有森林图 ────────────────────────────────────────────────
save_plot <- function(plot_obj, filename, w = 11, h = 6.5) {
  if (!is.null(plot_obj)) {
    ggsave(filename, plot = plot_obj, width = w, height = h,
           dpi = 300, bg = "white")
    message(sprintf("✓ Plot saved: %s", filename))
  }
}

# BBD plots
if (exists("res_bbd")) {
  save_plot(
    make_forest_plot(res_bbd,
                     "Truncating Variant Associations with BBD (Overall)",
                     "forest_BBD_overall.png"),
    "forest_BBD_overall.png"
  )
}
if (exists("res_nonprol")) {
  save_plot(
    make_forest_plot(res_nonprol,
                     "Truncating Variant Associations with Non-proliferative BBD",
                     "forest_BBD_nonproliferative.png"),
    "forest_BBD_nonproliferative.png"
  )
}
if (exists("res_prol")) {
  save_plot(
    make_forest_plot(res_prol,
                     "Truncating Variant Associations with Proliferative BBD (without atypia)",
                     "forest_BBD_proliferative.png"),
    "forest_BBD_proliferative.png"
  )
}

# DCIS / LCIS plots (if run)
if (exists("res_dcis")) {
  save_plot(
    make_forest_plot(res_dcis,
                     "Truncating Variant Associations with DCIS",
                     "forest_DCIS_truncating.png"),
    "forest_DCIS_truncating.png"
  )
}
if (exists("res_lcis")) {
  save_plot(
    make_forest_plot(res_lcis,
                     "Truncating Variant Associations with LCIS",
                     "forest_LCIS_truncating.png"),
    "forest_LCIS_truncating.png"
  )
}


# ── 4. 打印所有文件路径，方便你找到它们 ───────────────────────────────────────
message("\n", strrep("=", 55))
message("ALL EXPORTED FILES")
message(strrep("=", 55))

all_outputs <- c(
  xlsx_path,
  "forest_BBD_overall.png",
  "forest_BBD_nonproliferative.png",
  "forest_BBD_proliferative.png",
  "forest_DCIS_truncating.png",
  "forest_LCIS_truncating.png"
)

for (f in all_outputs) {
  full_path <- file.path(getwd(), f)
  exists_str <- if (file.exists(full_path)) "✓ exists" else "✗ not yet created"
  message(sprintf("  %-45s %s", f, exists_str))
}

message(sprintf("\nDirectory: %s", getwd()))
message("You can find all files in the Files tab (bottom-right panel in RStudio).")
message("Click 'More' > 'Show Folder in Finder/Explorer' to open the folder directly.")
