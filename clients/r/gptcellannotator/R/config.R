#' GPT Cell Annotator configuration
#'
#' @description
#' Creates or retrieves configuration used by `gptcellannotator` to talk to the
#' backend REST service or shell out to the CLI.
#'
#' @param base_url Base URL for the REST service. Defaults to the
#'   `GPTCA_BASE_URL` environment variable or `https://api.gpt-cell-annotator.org`.
#' @param api_key API key passed as `Authorization: Bearer` header when present.
#'   Defaults to `GPTCA_API_KEY`.
#' @param cli_path Optional path to the fallback CLI (`gca`). If `NULL`, the CLI
#'   is auto-detected using `Sys.which("gca")`.
#' @param timeout Timeout in seconds for HTTP requests.
#' @param retry_max Maximum retry attempts for transient HTTP failures.
#' @param retry_backoff Initial backoff in seconds for retries (exponential).
#' @param offline Logical flag forcing the CLI fallback.
#' @param user_agent Custom user agent string appended to requests.
#'
#' @return A `GptcaConfig` object.
#' @export
gptca_config <- function(
  base_url = Sys.getenv("GPTCA_BASE_URL", unset = "https://api.gpt-cell-annotator.org"),
  api_key = Sys.getenv("GPTCA_API_KEY", unset = NA_character_),
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  timeout = 120,
  retry_max = 3,
  retry_backoff = 1,
  offline = FALSE,
  user_agent = utils::packageName()
) {
  if (!nzchar(cli_path)) {
    cli_path <- Sys.which("gca")
    if (!nzchar(cli_path)) {
      cli_path <- NA_character_
    }
  }

  cfg <- list(
    base_url = normalize_base_url(base_url),
    api_key = if (!is.na(api_key) && nzchar(api_key)) api_key else NULL,
    cli_path = if (!is.null(cli_path) && nzchar(cli_path)) cli_path else NULL,
    timeout = validate_positive_number(timeout, "timeout"),
    retry_max = validate_non_negative_integer(retry_max, "retry_max"),
    retry_backoff = validate_positive_number(retry_backoff, "retry_backoff"),
    offline = isTRUE(offline),
    user_agent = gptca_user_agent(user_agent)
  )
  class(cfg) <- "GptcaConfig"
  cfg
}

#' @export
print.GptcaConfig <- function(x, ...) {
  cli::cli_text("{.strong GPT Cell Annotator configuration}")
  cli::cli_ul(c(
    "base_url" = x$base_url,
    "api_key" = ifelse(is.null(x$api_key), "<unset>", "<hidden>"),
    "cli_path" = ifelse(is.null(x$cli_path), "<not-found>", x$cli_path),
    "timeout" = sprintf("%.1f sec", x$timeout),
    "retry_max" = x$retry_max,
    "retry_backoff" = sprintf("%.1f sec", x$retry_backoff),
    "offline" = ifelse(isTRUE(x$offline), "yes", "no")
  ))
  invisible(x)
}

#' Activate a configuration
#'
#' @param config A `GptcaConfig` object.
#' @return The previously active configuration (invisibly).
#' @export
gptca_config_set <- function(config) {
  if (!inherits(config, "GptcaConfig")) {
    cli::cli_abort("{.arg config} must be a {.cls GptcaConfig}.")
  }
  prev <- gptca_config_get(default = NULL)
  assign("gptca_config", config, envir = .gptca_state)
  invisible(prev)
}

#' Retrieve the active configuration
#'
#' @param default Value returned if no config has been activated.
#' @return A `GptcaConfig` object.
#' @export
gptca_config_get <- function(default = gptca_config()) {
  cfg <- get0("gptca_config", envir = .gptca_state, ifnotfound = NULL)
  if (is.null(cfg)) {
    cfg <- default
    gptca_config_set(cfg)
  }
  cfg
}

#' Reset configuration to defaults
#'
#' @export
gptca_config_reset <- function() {
  assign("gptca_config", NULL, envir = .gptca_state)
  invisible(gptca_config_get())
}

.gptca_state <- new.env(parent = emptyenv())

normalize_base_url <- function(url) {
  if (!rlang::is_string(url) || !nzchar(url)) {
    cli::cli_abort("{.arg base_url} must be a non-empty string.")
  }
  url <- trimws(url)
  url <- sub("/+$", "", url)
  url
}

validate_positive_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value) || value <= 0) {
    cli::cli_abort("{.arg {arg}} must be a positive number.", arg = name)
  }
  as.numeric(value)
}

validate_non_negative_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0) {
    cli::cli_abort("{.arg {arg}} must be a non-negative number.", arg = name)
  }
  as.integer(value)
}

gptca_user_agent <- function(extra = NULL) {
  pkg <- tryCatch(
    utils::packageDescription("gptcellannotator"),
    error = function(...) NULL
  )
  version <- pkg[["Version"]] %||% "0.0.0"
  ua <- sprintf("gptcellannotator/%s", version)
  if (!is.null(extra) && nzchar(extra) && !identical(extra, utils::packageName())) {
    ua <- paste(ua, extra)
  }
  ua
}

`%||%` <- function(x, y) if (is.null(x)) y else x
