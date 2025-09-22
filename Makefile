.DEFAULT_GOAL := help

# default value if not provided
VERBOSE ?= 0
MODULE ?= *

.PHONY: help
help: ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*##"; printf "\033[36m\033[0m"} /^[a-zA-Z0-9_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)


check_for_sprocket:
	@echo "Checking if sprocket is available..."
	@if ! command -v sprocket >/dev/null 2>&1; then \
		echo >&2 "Error: sprocket is not installed or not in PATH. Install sprocket (https://sprocket.bio/installation.html)"; \
		exit 1; \
	else \
	  echo "sprocket version $$(sprocket --version | awk '{print $$2}')"; \
	fi;

check_for_uv:
	@echo "Checking if uv is available..."
	@if ! command -v uv >/dev/null 2>&1; then \
		echo >&2 "Error: uv is not installed or not in PATH. Install uv (https://docs.astral.sh/uv/getting-started/installation/)"; \
		exit 1; \
	else \
	  echo "uv version $$(uv --version | awk '{print $$2}')"; \
		uv run --python 3.13 --with miniwdl python -c "from importlib.metadata import version; print(f'miniwdl v{version(\"miniwdl\")}')"; \
	fi;

check_module:
	@if [ "$(MODULE)" != "*" ] && [ ! -d "modules/$(MODULE)" ]; then \
		echo >&2 "Error: Module '$(MODULE)' not found in modules/ directory"; \
		exit 1; \
	fi

##@ Linting

lint_sprocket: check_for_sprocket check_module ## Run sprocket lint on all modules or a specific module using MODULE=name
	@echo "Running sprocket lint..."
	@for dir in modules/$(MODULE)/; do \
		if [ -d "$$dir" ]; then \
			echo "Linting $$dir"; \
			sprocket lint "$$dir"; \
		fi; \
	done

lint_miniwdl: check_for_uv check_module ## Run miniwdl lint on all modules or a specific module using MODULE=name (use VERBOSE=1 for detailed output)
	@echo "Running miniwdl lint..."
	@for file in modules/$(MODULE)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Linting $$file"; \
			if [ "$(VERBOSE)" = "1" ]; then \
				uv run --python 3.13 --with miniwdl miniwdl check "$$file"; \
			else \
				uv run --python 3.13 --with miniwdl miniwdl check "$$file" >/dev/null; \
			fi; \
		fi; \
	done

lint: lint_sprocket lint_miniwdl ## Run all linting checks

##@ Run

run_sprocket: check_for_sprocket check_module ## Run sprocket run on all modules or a specific module using MODULE=name
	@echo "Running sprocket run..."
	@for file in modules/$(MODULE)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Testing sprocket run for $$file"; \
			entrypoint=$$(grep '^workflow' "$$file" | awk '{print $$2}'); \
			sprocket run "$$file" --entrypoint $$entrypoint; \
		fi; \
	done

run_miniwdl: check_for_uv check_module ## Run miniwdl run on all modules or a specific module using MODULE=name
	@echo "Running miniwdl run..."
	@for file in modules/$(MODULE)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "... for $$file"; \
			uv run --python 3.13 --with miniwdl miniwdl run "$$file"; \
		fi; \
	done

run: run_sprocket run_miniwdl ## Run all run checks
