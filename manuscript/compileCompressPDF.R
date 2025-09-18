setwd("manuscript")
knitr::knit2pdf("onStepsAndTurns.Rnw")
system("pdflatex Appendix_S1.tex")


system("gs -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -dNOPAUSE -dBATCH -sOutputFile=output.pdf onStepsAndTurns.pdf")
file.remove("onStepsAndTurns.pdf")
file.rename(from = "output.pdf", to = "onStepsAndTurns.pdf")

system("gs -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -dNOPAUSE -dBATCH -sOutputFile=output.pdf Appendix_S1.pdf")
file.remove("Appendix_S1.pdf")
file.rename(from = "output.pdf", to = "Appendix_S1.pdf")
setwd("..")
