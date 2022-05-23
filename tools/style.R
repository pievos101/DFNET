# Usage: Rscript tools/style.R
library(styler)
style <- function() {
    tidyverse_style(indent_by = 4)
}
style_dir("examples", style = style)
style_dir("R", style = style)
style_dir("tests", style = style)
