% Generated manually until roxygen can be run
\name{gptca_config}
\alias{gptca_config}
\alias{gptca_config_set}
\alias{gptca_config_get}
\alias{gptca_config_reset}
\title{Configure GPT Cell Annotator backend access}
\description{
Create or manage configuration objects used to communicate with the GPT Cell Annotator backend. The configuration stores base URLs, authentication, CLI fallbacks, and retry policies.
}
\usage{
gptca_config(
  base_url = Sys.getenv("GPTCA_BASE_URL", unset = "https://api.gpt-cell-annotator.org"),
  api_key = Sys.getenv("GPTCA_API_KEY", unset = NA_character_),
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  timeout = 120,
  retry_max = 3,
  retry_backoff = 1,
  offline = FALSE,
  user_agent = utils::packageName()
)

gptca_config_set(config)

gptca_config_get(default = gptca_config())

gptca_config_reset()
}
\arguments{
\item{base_url}{Base URL for the FastAPI backend.}

\item{api_key}{Optional API key passed as a Bearer token.}

\item{cli_path}{Path to the `gca` CLI for offline fallback; detected automatically when not supplied.}

\item{timeout}{Request timeout in seconds.}

\item{retry_max}{Maximum retry attempts for transient HTTP failures.}

\item{retry_backoff}{Initial delay in seconds applied to exponential backoff.}

\item{offline}{Logical flag forcing CLI mode.}

\item{user_agent}{Additional user-agent suffix appended to requests.}

\item{config}{A \verb{GptcaConfig} object to activate.}

\item{default}{Default configuration returned when none is active.}
}
\value{
\code{gptca_config()} returns an object of class \verb{GptcaConfig}. \code{gptca_config_set()} invisibly returns the previous config. \code{gptca_config_get()} returns the currently active config. \code{gptca_config_reset()} resets to defaults and returns the new config.
}
\examples{
if (interactive()) {
  cfg <- gptca_config(base_url = "https://api.gpt-cell-annotator.org", timeout = 60)
  gptca_config_set(cfg)
}
}
