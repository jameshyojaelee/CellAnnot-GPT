#' Internal: perform GPT Cell Annotator HTTP request
#'
#' @keywords internal
gptca_http_post <- function(path, body, config) {
  if (!inherits(config, "GptcaConfig")) {
    cli::cli_abort("{.arg config} must be a {.cls GptcaConfig}.")
  }

  url <- paste0(config$base_url, path)
  req <- httr2::request(url)
  req <- httr2::req_method(req, "POST")
  req <- httr2::req_body_json(req, data = body, auto_unbox = TRUE)
  req <- httr2::req_headers(
    req,
    Accept = "application/json",
    `Content-Type` = "application/json"
  )
  req <- httr2::req_user_agent(req, config$user_agent)
  req <- httr2::req_timeout(req, config$timeout)

  if (!is.null(config$api_key)) {
    req <- httr2::req_headers(req, Authorization = paste("Bearer", config$api_key))
  }

  req <- httr2::req_retry(
    req,
    max_tries = max(1L, config$retry_max + 1L),
    backoff = gptca_retry_backoff(config$retry_backoff),
    is_transient = gptca_is_transient_error
  )

  resp <- try(httr2::req_perform(req), silent = TRUE)
  if (inherits(resp, "try-error")) {
    cli::cli_abort(
      "Failed to contact GPT Cell Annotator service at {.uri {url}}.\n{conditionMessage(attr(resp, \"condition\"))}"
    )
  }
  httr2::resp_check_status(resp)
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

gptca_retry_backoff <- function(base_delay) {
  function(req) {
    tries <- req$policies$retry_tries %||% 1L
    delay <- base_delay * 2^(tries - 1L)
    stats::runif(1L, min = delay * 0.8, max = delay * 1.2)
  }
}

gptca_is_transient_error <- function(resp) {
  status <- httr2::resp_status(resp)
  status >= 500 || status == 429
}
