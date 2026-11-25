# Simulates the demo dataset in inst/extdata/. All sample IDs, genotypes, covariates and phenotypes are
# simulated. Run from the package root: Rscript data-raw/simulate_demo_data.R

set.seed(2026)
out <- "inst/extdata"
n <- 5000
m <- 500

## genotypes: independent SNPs, 0.5% missing calls
maf <- runif(m, 0.05, 0.5)
G <- sapply(maf, function(f) rbinom(n, 2, f))
G_obs <- G
G_obs[sample(length(G), 0.005 * length(G))] <- NA

## covariates
age <- round(runif(n, 40, 70))
sex <- rbinom(n, 1, 0.5)
pc1 <- rnorm(n)
pc2 <- rnorm(n)
cov <- 0.02 * (age - 55) + 0.2 * sex + 0.1 * pc1

## five causal SNPs; the phenotype is observed in 30% of the samples
causal <- c(50, 150, 250, 350, 450)
Gc <- scale(G[, causal])
labeled <- rbinom(n, 1, 0.3) == 1
observe <- function(y) ifelse(labeled, y, NA)

## continuous phenotype and its synthetic phenotype
y1 <- cov + c(Gc %*% c(0.11, -0.11, 0.10, 0.10, -0.10)) + rnorm(n)
s1 <- 0.2 + 0.8 * y1 + rnorm(n, 0, 0.35)

## binary phenotype (0/1); the synthetic phenotype is a predicted probability
y2 <- rbinom(n, 1, plogis(-1 + cov + c(Gc %*% c(0.24, -0.24, 0.22, 0.22, -0.22))))
s2 <- plogis(-1 + 2 * y2 + rnorm(n, 0, 1.2))

## count phenotype; the synthetic phenotype is a predicted count
y3 <- rpois(n, exp(1 + 0.5 * cov + c(Gc %*% c(0.07, -0.07, 0.06, 0.06, -0.06))))
s3 <- (y3 + 0.5) * exp(rnorm(n, 0, 0.3))

## PLINK files (A1 is the allele counted in G)
iid <- paste0("id", sprintf("%05d", sample(99999, n)))
fam <- data.frame(iid, iid, 0, 0, 0, -9)
bim <- data.frame(1, paste0("snp", seq_len(m)), 0, 1000 * seq_len(m), "A", "G")
write.table(fam, file.path(out, "demo.fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bim, file.path(out, "demo.bim"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

code <- c(3L, 2L, 0L)                       # A1 count 0, 1, 2 -> bed codes 11, 10, 00; missing -> 01
con <- file(file.path(out, "demo.bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
for (j in seq_len(m)) {
  b <- code[G_obs[, j] + 1L]
  b[is.na(b)] <- 1L
  b <- matrix(c(b, rep(0L, (4 - n %% 4) %% 4)), nrow = 4)
  writeBin(as.raw(colSums(b * c(1L, 4L, 16L, 64L))), con)
}
close(con)

## covariate and phenotype files
write.table(data.frame(IID = iid, age, sex, pc1 = round(pc1, 4), pc2 = round(pc2, 4)),
            file.path(out, "demo_covar.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
pheno <- function(y, s, name) {
  write.table(data.frame(IID = iid, pheno = observe(y), s_pheno = round(s, 4)),
              file.path(out, name), quote = FALSE, row.names = FALSE, sep = "\t")
}
pheno(round(y1, 4), s1, "demo_pheno_continuous.txt")
pheno(y2, s2, "demo_pheno_binary.txt")
pheno(y3, s3, "demo_pheno_count.txt")
