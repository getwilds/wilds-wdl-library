.DEFAULT_GOAL := help

# default value if not provided
VERBOSE ?= 0

.PHONY: help
help: ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*##"; printf "\033[36m\033[0m"} /^[a-zA-Z0-9_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Linting

lint_sprocket: ## Run sprocket lint on all WDL modules
	@echo "Running sprocket lint on all modules..."
	@echo "Checking if sprocket is available..."
	@if ! command -v sprocket >/dev/null 2>&1; then \
		echo >&2 "Error: sprocket is not installed or not in PATH. Install sprocket (https://sprocket.bio/installation.html)"; \
		exit 1; \
	fi
	@for dir in modules/*/; do \
		if [ -d "$$dir" ]; then \
			echo "Linting $$dir"; \
			sprocket lint "$$dir"; \
		fi; \
	done

lint_miniwdl: ## Run miniwdl lint on all WDL modules (use VERBOSE=1 for detailed output)
	@echo "Running miniwdl lint on all modules..."
	@echo "Checking if uv is available..."
	@if ! command -v uv >/dev/null 2>&1; then \
		echo >&2 "Error: uv is not installed or not in PATH. Install uv (https://docs.astral.sh/uv/getting-started/installation/)"; \
		exit 1; \
	fi
	@for file in modules/*/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Linting $$file"; \
			if [ "$(VERBOSE)" = "1" ]; then \
				uv run --python 3.13.5 --with miniwdl miniwdl check "$$file"; \
			else \
 	            result=`uv run --python 3.13.5 --with miniwdl miniwdl check "$$file"` || echo $$result; \
     	fi; \
		fi; \
	done

lint: lint_sprocket lint_miniwdl ## Run all linting checks
