`htmlresult` <-
function(hgOver, filename, app = FALSE, categorySize = 1, pvalue = 0.05) {
  write(paste("<br>------<br>", length(geneIdUniverse(hgOver, cond = conditional(hgOver))), testName(hgOver)[1], testName(hgOver)[2], "ids tested(", dim(summary(hgOver, pvalue = 0.05))[1], " p<0.05)<br>"), file = filename, append = app)
  write(paste("Gene universe size:", universeMappedCount(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Selected gene set size:", geneMappedCount(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Conditional:", conditional(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Annotation:", annotation(hgOver), "<br><br>"), file = filename, append = TRUE)
  htmlReport(hgOver, summary.args = list(htmlLinks = TRUE, categorySize = categorySize, pvalue = pvalue), file = filename, append = TRUE)
}

