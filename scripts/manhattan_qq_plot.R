# ============================================
# Manhattan Plot and QQ Plot
# GWAS results visualization
# Author: Azadeh Golduzian | UNM PhD Research
# ============================================
R --vanilla -q <<'EOF'
library(data.table)
library(ggplot2)

BASE <- "/home/jupyter/Chapter2_Azadeh_01072026"
CODE <- file.path(BASE, "Code_produced")
FIG  <- file.path(BASE, "Figures")
BED  <- file.path(BASE, "TRANS_windows_0.5Mb.nochr.bed")

# ----------------------------
# Load TRANS regions
# ----------------------------
stopifnot(file.exists(BED))
tr <- fread(BED, header=FALSE)
setnames(tr, c("CHR","TSTART","TEND"))
tr[, CHR := as.integer(CHR)]
setkey(tr, CHR, TSTART, TEND)

# ----------------------------
# Helpers
# ----------------------------
read_gwas <- function(f, label){
  if(!file.exists(f)) stop("Missing GWAS file: ", f)
  dt <- fread(f)
  if("#CHROM" %in% names(dt)) setnames(dt, "#CHROM", "CHR")
  if(!"CHR" %in% names(dt)) setnames(dt, 1, "CHR")

  dt <- dt[TEST=="ADD" & !is.na(P) & P > 0]
  dt[, CHR := as.integer(gsub("^chr","", as.character(CHR)))]
  dt <- dt[CHR >= 1L & CHR <= 22L]
  dt[, .(CHR, BP=as.numeric(POS), ID, P, neglog10p=-log10(P), group=label)]
}

calc_lambda <- function(p){
  # genomic inflation factor (lambda GC)
  median(qchisq(1 - p, df=1), na.rm=TRUE) / qchisq(0.5, df=1)
}

tag_reg <- function(g){
  g[, rid := .I]
  g[, `:=`(QS=BP, QE=BP)]
  setkey(g, CHR, QS, QE)

  ov <- foverlaps(g, tr,
                  by.x=c("CHR","QS","QE"),
                  by.y=c("CHR","TSTART","TEND"),
                  nomatch=NA)

  reg_ids <- ov[!is.na(TSTART), unique(rid)]
  g[, region := fifelse(rid %in% reg_ids, "REG", "OTHER")]
  g[, c("rid","QS","QE") := NULL]
  g[]
}

add_cumpos <- function(g){
  chrlen <- g[, .(maxBP=max(BP, na.rm=TRUE)), by=CHR][order(CHR)]
  chrlen[, offset := cumsum(as.numeric(shift(maxBP, fill=0)))]
  g <- merge(g, chrlen[,.(CHR, offset)], by="CHR", all.x=TRUE)
  g[, pos := BP + offset]
  centers <- chrlen[, .(CHR, center=offset + maxBP/2)]
  list(g=g, centers=centers)
}

# ----------------------------
# Font: try Arial, else fallback to sans
# ----------------------------
pick_font <- function(){
  # pdfFonts() often includes "Helvetica" etc; Arial may not exist on Linux.
  fams <- names(grDevices::pdfFonts())
  if ("Arial" %in% fams) return("Arial")
  return("sans")
}
FONT_FAMILY <- pick_font()

# ----------------------------
# Miami plot (grant style)
# ----------------------------
make_miami_grant <- function(g1, g2, label1, label2, title, outfile_pdf, outfile_png){

  g1[, group := label1]
  g2[, group := label2]
  g <- rbindlist(list(g1, g2), use.names=TRUE, fill=TRUE)

  g <- tag_reg(g)

  tmp <- add_cumpos(g)
  g <- tmp$g
  centers <- tmp$centers

  # show fewer chromosome labels to avoid overlap
  show_chr <- c(1:10, 12, 14, 16, 18, 20, 22)
  centers_show <- centers[CHR %in% show_chr]

  lam1 <- calc_lambda(g[group==label1, P])
  lam2 <- calc_lambda(g[group==label2, P])
  n_sig1 <- sum(g[group==label1, P < 5e-8], na.rm=TRUE)
  n_sig2 <- sum(g[group==label2, P < 5e-8], na.rm=TRUE)

  g[, is_sig := P < 5e-8]
  g[, y_miami := fifelse(group==label1,  neglog10p, -neglog10p)]

  # Extend y-axis range so points don't hit the edge (Davorka request)
  thr  <- -log10(5e-8)
  ymax <- max(abs(g$y_miami), na.rm=TRUE) * 1.15   # 15% padding
  ylim <- c(-ymax, ymax)

  p <- ggplot(g, aes(pos, y_miami)) +
    geom_point(
      data=g[is_sig==FALSE],
      aes(color=group, shape=region),
      size=0.8, alpha=0.30
    ) +
    geom_point(
      data=g[is_sig==TRUE],
      aes(color=group, shape=region),
      size=2.5, alpha=0.90
    ) +
    geom_hline(yintercept=0, linewidth=0.4, color="black") +
    geom_hline(yintercept=c(thr, -thr), linetype="dashed", linewidth=0.3, color="gray40") +
    scale_x_continuous(breaks=centers_show$center, labels=centers_show$CHR, expand=c(0.01,0)) +
    scale_y_continuous(limits=ylim, expand=c(0.02,0)) +
    scale_shape_manual(values=c("OTHER"=16, "REG"=17), labels=c("Other", "Regulatory")) +
    scale_color_manual(values=setNames(c("#0072B2", "#D55E00"), c(label1, label2))) +
    labs(
      x="Chromosome",
      y=expression(paste(-log[10], "(", italic(P), ")")),
      title=title,
      subtitle=paste0("λ: ", label1, "=", round(lam1,2),
                      ", ", label2, "=", round(lam2,2),
                      " | GW-sig: ", label1, "=", n_sig1,
                      ", ", label2, "=", n_sig2),
      shape="Region",
      color="Group"
    ) +
    theme_classic(base_size=9, base_family=FONT_FAMILY) +
    theme(
      legend.position="top",
      legend.box="horizontal",
      legend.margin=margin(0,0,0,0),
      legend.key.size=unit(0.30, "cm"),
      legend.text=element_text(size=7),
      legend.title=element_text(size=7),
      plot.title=element_text(face="bold", hjust=0.5, size=9),
      plot.subtitle=element_text(hjust=0.5, size=6),
      axis.text.x=element_text(size=6),
      axis.text.y=element_text(size=7),
      axis.title=element_text(size=8),
      plot.margin=margin(2,5,2,2)
    ) +
    guides(
      color=guide_legend(order=1, override.aes=list(size=2)),
      shape=guide_legend(order=2, override.aes=list(size=2))
    )

  # ---- PDF via cairo (Davorka: cairo better than ggsave for PDF) ----
  grDevices::cairo_pdf(outfile_pdf, width=4, height=2.5)
  print(p)
  grDevices::dev.off()
  cat("Saved PDF:", outfile_pdf, "\n")

  # ---- PNG 600 dpi ----
  ggsave(outfile_png, p, width=4, height=2.5, dpi=600, device="png")
  cat("Saved PNG:", outfile_png, "\n")
  cat("Font family used:", FONT_FAMILY, "\n")
}

# ==============================================================
# 1) DBP: Young vs Old
# ==============================================================
cat("\n=== Panel 1: DBP Young vs Old ===\n")
young_dbp <- read_gwas(file.path(CODE, "Under60", "gwas_u60_dbp_INT_20PC.dbp_INT.glm.linear"), "Young (<60)")
old_dbp   <- read_gwas(file.path(CODE, "Over60",  "gwas_o60_dbp_INT_20PC.dbp_INT.glm.linear"), "Old (≥60)")

make_miami_grant(
  young_dbp, old_dbp,
  "Young (<60)", "Old (≥60)",
  "DBP: Young vs Old",
  file.path(FIG, "dbp_miami_age.pdf"),
  file.path(FIG, "dbp_miami_age.png")
)

# ==============================================================
# 2) DBP Young: Male vs Female
# ==============================================================
cat("\n=== Panel 2: DBP Young - Male vs Female ===\n")
young_m_dbp <- read_gwas(file.path(CODE, "Under60", "gwas_u60_M_dbp_INT_5PC.dbp_INT.glm.linear"), "Male")
young_f_dbp <- read_gwas(file.path(CODE, "Under60", "gwas_u60_F_dbp_INT_5PC.dbp_INT.glm.linear"), "Female")

make_miami_grant(
  young_m_dbp, young_f_dbp,
  "Male", "Female",
  "DBP Young (<60): Male vs Female",
  file.path(FIG, "dbp_young_sex.pdf"),
  file.path(FIG, "dbp_young_sex.png")
)

# ==============================================================
# 3) DBP Old: Male vs Female
# ==============================================================
cat("\n=== Panel 3: DBP Old - Male vs Female ===\n")
old_m_dbp <- read_gwas(file.path(CODE, "Over60", "gwas_o60_M_dbp_INT_5PC.dbp_INT.glm.linear"), "Male")
old_f_dbp <- read_gwas(file.path(CODE, "Over60", "gwas_o60_F_dbp_INT_5PC.dbp_INT.glm.linear"), "Female")

make_miami_grant(
  old_m_dbp, old_f_dbp,
  "Male", "Female",
  "DBP Old (≥60): Male vs Female",
  file.path(FIG, "dbp_old_sex.pdf"),
  file.path(FIG, "dbp_old_sex.png")
)

cat("\nDone! Saved 3 panels.\n")
cat("Target size: 4 x 2.5 inches\n")
cat("PNG resolution: 600 dpi\n")
cat("PDF: cairo_pdf\n")
cat("Requested font: Arial 9 (fallback if not installed)\n")
EOF
