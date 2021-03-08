assignments for experimental design and data analysis

some notes: 

1. use `rmarkdown::render("*.Rmd")` to render rmd files
2. use `<-` for value assignment
3. use `with()` instead of `attach()` to avoid naming conflicts since we working on one rmd file
    - similarly, use explicitly qualify like `car::avPlots()` instead of `library(car)`, except that we use the library a lot throughtout the assignment
4. use `fig.height` and `fig.width` to control the size of the figure, as a rule of thumb, 2.5 is good for `fig.height` of one figure
    - so, `{r, fig.height=2.5, fig.width=6, fig.align='center'}` for two figures in the same row
5. use `mar` parameter in `par()` to shrink the blank spaces, e.g. `par(mfrow=c(1, 2), mar=c(3, 5, 1, 0))`
6. break lines to make it easier for reviewing one certain sentence
