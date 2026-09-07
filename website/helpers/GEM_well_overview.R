escape_html <- function(value) {
  htmltools::htmlEscape(as.character(value), attribute = FALSE)
}

format_GEM_well_value <- function(value, column) {
  if (is.na(value)) {
    return("<code class=\"gem-well-missing\">NA</code>")
  }

  if (identical(column, "GEM_well_QC_exclude_list")) {
    value <- gsub(" ;; ", " ;;<wbr> ", escape_html(value), fixed = TRUE)
    return(paste0("<code class=\"gem-well-QC-value\">", value, "</code>"))
  }

  paste0("<code>", escape_html(value), "</code>")
}

emit_GEM_well_demo_table <- function(
  GEM_well_config_file,
  dictionary_file
) {
  demo <- list(GEM_wells = utils::read.delim(
    GEM_well_config_file, check.names = FALSE, colClasses = "character"
  ))
  dictionary <- utils::read.delim(
    dictionary_file,
    check.names = FALSE,
    colClasses = "character",
    na.strings = character()
  )
  required_dictionary_columns <- c(
    "column",
    "label",
    "category",
    "purpose",
    "details",
    "is_required"
  )

  if (!identical(names(dictionary), required_dictionary_columns)) {
    stop(
      "GEM-well dictionary columns must be: ",
      paste(required_dictionary_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(dictionary$column)) {
    stop("GEM-well dictionary contains duplicate column names", call. = FALSE)
  }
  if (!all(dictionary$is_required %in% c("TRUE", "FALSE"))) {
    stop("GEM-well dictionary is_required values must be TRUE or FALSE", call. = FALSE)
  }
  if (!setequal(dictionary$column, names(demo$GEM_wells))) {
    missing_columns <- setdiff(names(demo$GEM_wells), dictionary$column)
    extra_columns <- setdiff(dictionary$column, names(demo$GEM_wells))
    stop(
      "GEM-well dictionary does not match cfg_GEM_wells.tsv. Missing: ",
      paste(missing_columns, collapse = ", "),
      "; extra: ",
      paste(extra_columns, collapse = ", "),
      call. = FALSE
    )
  }

  dictionary <- dictionary[match(names(demo$GEM_wells), dictionary$column), , drop = FALSE]
  header_cells <- vapply(
    seq_len(nrow(dictionary)),
    function(index) {
      field <- dictionary[index, ]
      tooltip_id <- paste0("gem-well-help-", gsub("[^[:alnum:]-]", "-", field$column))
      required_class <- if (identical(field$is_required, "TRUE")) {
        " class=\"tsv-required-column\""
      } else {
        ""
      }
      paste0(
        "<th scope=\"col\"",
        required_class,
        "><span class=\"gem-well-column-heading\"><code>",
        escape_html(field$column),
        "</code><span class=\"gem-well-tooltip\"><button type=\"button\" aria-describedby=\"",
        tooltip_id,
        "\" aria-label=\"More information about ",
        escape_html(field$label),
        "\">i</button><span id=\"",
        tooltip_id,
        "\" role=\"tooltip\"><strong>",
        escape_html(field$label),
        ".</strong> ",
        escape_html(field$purpose),
        " ",
        escape_html(field$details),
        "</span></span></span></th>"
      )
    },
    character(1)
  )
  body_rows <- vapply(
    seq_len(nrow(demo$GEM_wells)),
    function(row_index) {
      values <- vapply(
        names(demo$GEM_wells),
        function(column) {
          format_GEM_well_value(demo$GEM_wells[[column]][row_index], column)
        },
        character(1)
      )
      cell_classes <- ifelse(
        dictionary$is_required == "TRUE",
        " class=\"tsv-required-column\"",
        ""
      )
      paste0(
        "<tr>",
        paste0("<td", cell_classes, ">", values, "</td>", collapse = ""),
        "</tr>"
      )
    },
    character(1)
  )

  cat(
    "<div class=\"gem-well-overview\"><table>",
    "<caption>cfg_GEM_wells.tsv</caption>",
    "<thead><tr>",
    paste(header_cells, collapse = ""),
    "</tr></thead><tbody>",
    paste(body_rows, collapse = ""),
    "</tbody></table></div>\n",
    sep = ""
  )
}

