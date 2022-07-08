## DFNET_utils - Various utilities.

## Copyright © 2021 Bastian Pfeifer <bastianxpfeifer@gmail.com>
## Copyright © 2022 Liliana Prikler <liliana.prikler@ist.tugraz.at>

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

#' Extract the multi-modal target vector from \code{features}
#'
#' @param features potentially multi-modal features as a matrix
#' @return A list of shape \code{(target, is.multi_modal)}, where \code{target}
#' is the target vector and \code{is.multi_modal} is true if the features are
#' multi-modal.
multi_modal_target <- function(features) {
    is.multi_modal <- !is.data.frame(features)
    if (is.multi_modal) {
        target <- features[[1]][, "target"]
    } else {
        target <- features[, "target"]
    }
    return(list(target = target, is.multi_modal = is.multi_modal))
}
