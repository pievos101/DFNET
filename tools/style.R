# Usage: Rscript tools/style.R
style <- function() {
    styler::tidyverse_style(indent_by = 4)
}
styler::style_dir("examples", style = style)
styler::style_dir("R", style = style)
styler::style_dir("tests", style = style)
