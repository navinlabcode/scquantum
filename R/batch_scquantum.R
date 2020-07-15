#' Batch execution of scQuantum single-cell ploidy calculation
#'
#' This function runs scQuantum ploidy inference in batches of cells. Allows parallel usage with the help of \code{future}.
#' @param df A data frame containing bincounts from DNA sequenced single-cells.
#' @param chrom A vector containing chromosome information for the bincounts.
#' @param start A vector containing start position information for the bincounts.
#' @param end A vector containing end position information for the bincounts.
#' @param penalty Numeric passed on to \link{\code{scquantum::ploidy.inference}}. Defaults to scquantum default value of 25.
#' @import purrr
#' @importFrom furrr future_map
#' @examples
#'
#'
#' @return A list object containing the scQuantum results for each cell.
#' @export
batch_scquantum <- function(df,
                            chrom = NULL,
                            start = NULL,
                            end = NULL,
                            penalty = 25) {
  furrr::future_map(df,
                    scquantum::ploidy.inference,
                    chrom,
                    start,
                    end,
                    penalty = penalty)

}

#' Accessor function to obtain the ploidy values of a scQuantum batch list
#'
#' This function runs scQuantum ploidy inference in batches of cells. Allows parallel usage with the help of \code{future}.
#' @param scquantum_list A list object containing the results from batch ploidy
#' @import purrr
#' @importFrom purrr map
#' @importFrom purrr pluck
#' @importFrom dplyr bind_rows
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename
#' @examples
#'
#'
#' @return A a data frame with the scquantum infered ploidy of each cell.
#' @export
batch_ploidy <- function(scquantum_list) {
  ploidies_list <- purrr::map(scquantum_list,
                              purrr::pluck,
                              "ploidy")

  ploidies_df <- dplyr::bind_rows(ploidies_list) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(ploidy = "V1")

  return(ploidies_df)

}

#' Accessor function to obtain the confidence values of a scQuantum batch list
#'
#' This function runs scQuantum confidence of the ploidy inference in batches of cells. Allows parallel usage with the help of \code{future}.
#' @param scquantum_list A list object containing the results from batch ploidy
#' @import purrr

#' @examples
#'
#'
#' @return A a data frame with the scquantum peak heights of each cell.
#' @export
batch_confidence <- function(scquantum_list) {
  confidence_list <- purrr::map(scquantum_list,
                                purrr::pluck,
                                "peak_height")

  confidence_df <- dplyr::bind_rows(ploidies_list) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(peak_height = "V1")

  return(ploidies_df)

}
