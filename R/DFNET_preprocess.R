## DFNET_preprocess - Preprocess input graphs.

## Copyright © 2021, 2022 Bastian Pfeifer <bastianxpfeifer@gmail.com>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

DFNET_preprocess <- function(DFNET_graph) {
    # XXX: Express in terms of relat
    range01 <- function(x) {
        res <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        if (is.infinite(res[1]) || is.na(res[1])) {
            return(x)
        } else {
            return(res)
        }
    }


    FEAT <- DFNET_graph$feature.matrix

    for (xx in 1:length(FEAT)) {
        FEAT[[xx]] <- as.data.frame(apply(FEAT[[xx]], 2, range01))
        colnames(FEAT[[xx]]) <- colnames(DFNET_graph$feature.matrix[[xx]])
    }

    DFNET_graph$feature.matrix <- FEAT
    return(DFNET_graph)
}
