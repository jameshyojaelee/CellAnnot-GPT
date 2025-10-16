# Simple helpers for gptcellannotator R package

PKG_PATH := gptcellannotator

.PHONY: install test check lint format doc pkgdown clean

install:
	Rscript -e "devtools::install('$(PKG_PATH)', upgrade = 'never')"

test:
	Rscript -e "devtools::test('$(PKG_PATH)')"

check:
	Rscript -e \"devtools::check('$(PKG_PATH)', error_on='warning', args=c('--as-cran'))\"

lint:
	Rscript -e "lintr::lint_package('$(PKG_PATH)')"

format:
	Rscript -e "styler::style_pkg('$(PKG_PATH)')"

doc:
	Rscript -e "devtools::document('$(PKG_PATH)')"

pkgdown:
	Rscript -e \"pkgdown::build_site('$(PKG_PATH)')\"

clean:
	Rscript -e "devtools::clean_dll('$(PKG_PATH)'); unlink(list.files('$(PKG_PATH)', pattern='\\\\.Rcheck$', full.names = TRUE), recursive = TRUE)"
