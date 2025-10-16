test_that("gptca_config uses env defaults", {
  withr::with_envvar(
    c("GPTCA_BASE_URL" = "https://example.org/api", "GPTCA_API_KEY" = "secret"),
    {
      cfg <- gptca_config()
      expect_s3_class(cfg, "GptcaConfig")
      expect_equal(cfg$base_url, "https://example.org/api")
      expect_equal(cfg$api_key, "secret")
    }
  )
})

test_that("gptca_config_set stores config", {
  cfg <- gptca_config(base_url = "https://example.org")
  prev <- gptca_config_set(cfg)
  expect_s3_class(prev, "GptcaConfig")
  expect_identical(gptca_config_get(), cfg)
})
