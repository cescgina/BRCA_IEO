```{r volcano3, echo=FALSE, out.width="600px", fig.cap="Figure S16: Volcano plot of the results of the DE analysis for the model with paired design."}
par(mar = c(4, 5, 2, 2), cex = 1, mfrow = c(1,2))
volcanoplot(fit, coef = 2, highlight = 7, fit$genes$symbol, main = "Model 1", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
volcanoplot(fitpaired, coef = 2, highlight = 7, fitpaired$genes$symbol, main = "Model 2", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
```
