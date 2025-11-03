.DEFAULT_GOAL := help

# default values if not provided
VERBOSE ?= 0
NAME ?= *
WOMTOOL ?= 86
WOMTOOL_JAR ?= womtool-$(WOMTOOL).jar
CROMWELL ?= 86
CROMWELL_JAR ?= cromwell-$(CROMWELL).jar

.PHONY: help
help: ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*##"; printf "\033[36m\033[0m"} /^[a-zA-Z0-9_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)


check_sprocket:
	@echo "Checking if sprocket is available..."
	@if ! command -v sprocket >/dev/null 2>&1; then \
		echo >&2 "Error: sprocket is not installed or not in PATH. Install sprocket (https://sprocket.bio/installation.html)"; \
		exit 1; \
	else \
	  echo "sprocket version $$(sprocket --version | awk '{print $$2}')"; \
	fi;

check_uv:
	@echo "Checking if uv is available..."
	@if ! command -v uv >/dev/null 2>&1; then \
		echo >&2 "Error: uv is not installed or not in PATH. Install uv (https://docs.astral.sh/uv/getting-started/installation/)"; \
		exit 1; \
	else \
	  echo "uv version $$(uv --version | awk '{print $$2}')"; \
		uv run --python 3.13 --with miniwdl python -c "from importlib.metadata import version; print(f'miniwdl v{version(\"miniwdl\")}')"; \
	fi;

check_wdlparse:
	@echo "Checking if wdlparse is available..."
	@if ! command -v wdlparse >/dev/null 2>&1; then \
		echo >&2 "Error: wdlparse is not installed or not in PATH. Install wdlparse (https://github.com/getwilds/wdlparse?tab=readme-ov-file#from-releases)"; \
		exit 1; \
	else \
	  echo "wdlparse version $$(wdlparse --version | awk '{print $$2}')"; \
	fi;

check_name:
	@if [ "$(NAME)" != "*" ]; then \
		if [ ! -d "modules/$(NAME)" ] && [ ! -d "vignettes/$(NAME)" ] && [ ! -d "workflows/$(NAME)" ]; then \
			echo >&2 "Error: '$(NAME)' not found in modules/, vignettes/, or workflows/ directory"; \
			exit 1; \
		fi; \
	fi

check_java:
	@echo "Checking your java version..."
	@if ! command -v java >/dev/null 2>&1; then \
		echo >&2 "Error: java is not installed or not in PATH."; \
		echo >&2 "Install Java 17 or 21 (showing Java 21 instructions as an example):"; \
		echo >&2 "  - macOS: brew install openjdk@21 && sudo ln -sfn /opt/homebrew/opt/openjdk@21/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-21.jdk"; \
		echo >&2 "  - Ubuntu/Debian: sudo apt install openjdk-21-jdk"; \
		echo >&2 "  - Other: https://adoptium.net/?variant=openjdk21&jvmVariant=hotspot"; \
		exit 1; \
	else \
		java_version=$$(java -version 2>&1 | head -n 1 | awk -F '"' '{print $$2}' | awk -F '.' '{print $$1}'); \
		if [ "$$java_version" = "17" ] || [ "$$java_version" = "21" ]; then \
			echo "Java $$java_version is installed and compatible with Cromwell"; \
		else \
			echo >&2 "Error: Default Java version is $$java_version, but Java 17 or 21 is required for Cromwell."; \
			echo >&2 "Install Java 17 or 21 and set it as default (showing Java 21 instructions as an example):"; \
			echo >&2 "  - macOS: brew install openjdk@21 && sudo ln -sfn /opt/homebrew/opt/openjdk@21/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-21.jdk"; \
			echo >&2 "  - Ubuntu/Debian: sudo apt install openjdk-21-jdk && sudo update-alternatives --set java /usr/lib/jvm/java-21-openjdk-amd64/bin/java"; \
			echo >&2 "  - Other: https://adoptium.net/?variant=openjdk21&jvmVariant=hotspot"; \
			exit 1; \
		fi; \
	fi;

check_womtool:
	@echo "Checking if WOMtool is available..."
	@if [ ! -f "$(WOMTOOL_JAR)" ]; then \
		echo "... $(WOMTOOL_JAR) not found, attempting to download..."; \
		if wget -O "$(WOMTOOL_JAR)" -q https://github.com/broadinstitute/cromwell/releases/download/$(WOMTOOL)/$(WOMTOOL_JAR); then \
			echo "... $(WOMTOOL_JAR) downloaded successfully"; \
		else \
			echo >&2 "... Error: Failed to download $(WOMTOOL_JAR). Please download it manually from https://github.com/broadinstitute/cromwell/releases"; \
			exit 1; \
		fi; \
	else \
		echo "... $(WOMTOOL_JAR) found"; \
	fi;

check_cromwell:
	@echo "Checking if cromwell is available..."
	@if [ ! -f "$(CROMWELL_JAR)" ]; then \
		echo "... $(CROMWELL_JAR) not found, attempting to download..."; \
		if wget -O "$(CROMWELL_JAR)" -q https://github.com/broadinstitute/cromwell/releases/download/$(CROMWELL)/$(CROMWELL_JAR); then \
			echo "... $(CROMWELL_JAR) downloaded successfully"; \
		else \
			echo >&2 "... Error: Failed to download $(CROMWELL_JAR). Please download it manually from https://github.com/broadinstitute/cromwell/releases"; \
			exit 1; \
		fi; \
	else \
		echo "... $(CROMWELL_JAR) found"; \
	fi;

##@ Linting

lint_sprocket: check_sprocket check_name ## Run sprocket lint on modules, vignettes, and workflows (use NAME=foo for specific item)
	@echo "Running sprocket lint..."
	@for dir in modules/$(NAME)/ vignettes/$(NAME)/ workflows/$(NAME)/; do \
		if [ -d "$$dir" ]; then \
			echo "Linting $$dir"; \
			sprocket lint "$$dir"; \
		fi; \
	done

lint_miniwdl: check_uv check_name ## Run miniwdl lint on modules, vignettes, and workflows (use NAME=foo, VERBOSE=1 for details)
	@echo "Running miniwdl lint..."
	@for file in modules/$(NAME)/*.wdl vignettes/$(NAME)/*.wdl workflows/$(NAME)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Linting $$file"; \
			if [ "$(VERBOSE)" = "1" ]; then \
				uv run --python 3.13 --with miniwdl miniwdl check "$$file"; \
			else \
				uv run --python 3.13 --with miniwdl miniwdl check "$$file" >/dev/null; \
			fi; \
		fi; \
	done

lint_womtool: check_java check_womtool check_name ## Run WOMtool validate on modules, vignettes, and workflows (use NAME=foo for specific item)
	@echo "Running WOMtool validate..."
	@set -e; for file in modules/$(NAME)/*.wdl vignettes/$(NAME)/*.wdl workflows/$(NAME)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Validating $$file"; \
			java -jar $(WOMTOOL_JAR) validate "$$file"; \
		fi; \
	done


lint: lint_sprocket lint_miniwdl lint_womtool ## Run all linting checks

##@ Run

run_sprocket: check_sprocket check_name check_wdlparse ## Run sprocket on testrun.wdl files (use NAME=foo for specific item)
	@echo "Running sprocket on testrun.wdl files..."
	@set -e; for dir in modules/$(NAME)/ vignettes/$(NAME)/ workflows/$(NAME)/; do \
		if [ -d "$$dir" ] && [ -f "$$dir/testrun.wdl" ]; then \
			echo "... Running $$(basename $$dir)"; \
			entrypoint=$$(wdlparse parse --format json "$$dir/testrun.wdl" | jq -r '.wdl.workflows[].name'); \
			echo "... Using entrypoint: $$entrypoint"; \
			sprocket run "$$dir/testrun.wdl" --entrypoint $$entrypoint; \
		fi; \
	done

run_miniwdl: check_uv check_name ## Run miniwdl on testrun.wdl files (use NAME=foo for specific item)
	@echo "Running miniwdl on testrun.wdl files..."
	@set -e; for dir in modules/$(NAME)/ vignettes/$(NAME)/ workflows/$(NAME)/; do \
		if [ -d "$$dir" ] && [ -f "$$dir/testrun.wdl" ]; then \
			echo "... Running $$(basename $$dir)"; \
			uv run --python 3.13 --with miniwdl miniwdl run "$$dir/testrun.wdl"; \
		fi; \
	done

run_cromwell: check_java check_cromwell check_name ## Run Cromwell on testrun.wdl files (use NAME=foo for specific item)
	@echo "Running Cromwell on testrun.wdl files..."
	@set -e; for dir in modules/$(NAME)/ vignettes/$(NAME)/ workflows/$(NAME)/; do \
		if [ -d "$$dir" ] && [ -f "$$dir/testrun.wdl" ]; then \
			echo "... Running $$(basename $$dir)"; \
			java -jar $(CROMWELL_JAR) run "$$dir/testrun.wdl"; \
		fi; \
	done

run: run_sprocket run_miniwdl run_cromwell ## Run all engines on testrun.wdl files

##@ Documentation

docs-preview: check_sprocket check_uv ## Build and serve documentation preview locally
	@echo "Building documentation preview..."
	@echo "Step 1/7: Checking for uncommitted changes..."
	@if git diff --quiet HEAD -- modules/ vignettes/ workflows/ && \
	   ! git ls-files --others --exclude-standard modules/ vignettes/ workflows/ | grep -q .; then \
		echo "... No uncommitted changes detected"; \
	else \
		echo ""; \
		echo "WARNING: You have uncommitted changes in modules/, vignettes/, or workflows/"; \
		echo "         The documentation will be built from your LAST COMMITTED version."; \
		echo "         Your uncommitted changes will be preserved but NOT included in the preview."; \
		echo "         To include them, commit your changes first, then run 'make docs-preview'."; \
		echo ""; \
		read -p "Continue anyway? [y/N] " -n 1 -r; \
		echo ""; \
		if [[ ! $$REPLY =~ ^[Yy]$$ ]]; then \
			echo "Aborted."; \
			exit 1; \
		fi; \
	fi
	@echo "Step 2/7: Saving current state..."
	@git stash push -u -m "docs-preview: temporary stash before building docs" -- modules/ vignettes/ workflows/
	@echo "Step 3/7: Creating preambles..."
	@uv run --python 3.13 .github/scripts/make_preambles.py
	@echo "Step 4/7: Creating .sprocketignore..."
	@printf '%s\n' '# Excluding test run WDL'\''s from documentation builds' 'modules/**/testrun.wdl' 'vignettes/**/testrun.wdl' 'workflows/**/testrun.wdl' > .sprocketignore
	@echo "Step 5/7: Building docs with sprocket..."
	@sprocket dev doc -v --homepage README.md --logo WILDSWDLNameLogo.svg .
	@echo "Step 6/7: Post-processing documentation..."
	@uv run --python 3.13 .github/scripts/postprocess_docs.py
	@echo "Step 7/7: Restoring original state..."
	@git restore modules/ vignettes/ workflows/ 2>/dev/null || true
	@rm -f .sprocketignore
	@if git stash list | grep -q "docs-preview: temporary stash before building docs"; then \
		echo "... Restoring your previous changes"; \
		git stash pop --index 2>/dev/null || git stash pop 2>/dev/null || true; \
	fi
	@echo ""
	@echo "Documentation built successfully!"
	@echo "To preview, run: make docs-serve"
	@echo "Or open: docs/index.html"

docs-serve: ## Serve documentation locally (requires docs to be built first)
	@echo "Starting local documentation server..."
	@if [ ! -d "docs" ]; then \
		echo >&2 "Error: docs directory not found. Run 'make docs-preview' first."; \
		exit 1; \
	fi
	@echo "Documentation will be available at: http://localhost:8000"
	@echo "Press Ctrl+C to stop the server"
	@cd docs && python3 -m http.server 8000

docs: docs-preview docs-serve ## Build and serve documentation in one command
