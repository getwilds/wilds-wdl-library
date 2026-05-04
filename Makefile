.DEFAULT_GOAL := help

# default values if not provided
VERBOSE ?= 0
NAME ?= *
WOMTOOL ?= 92
WOMTOOL_JAR ?= womtool-$(WOMTOOL).jar
CROMWELL ?= 92
CROMWELL_JAR ?= cromwell-$(CROMWELL).jar
MINIWDL ?= 1.13.0
SPROCKET_MIN ?= 0.22.0
SPROCKET_CONFIG ?=
SPROCKET_CONFIG_FLAG := $(if $(SPROCKET_CONFIG),-c $(SPROCKET_CONFIG),)
TYPE ?= all
NUM_RETRIES ?= 0
# TARGET=ci uses testrun.wdl; TARGET=hpc prefers testrun_hpc.wdl when present
# and falls back to testrun.wdl. Items lacking testrun.wdl are HPC-only and
# are skipped when TARGET=ci.
TARGET ?= ci

.PHONY: help
help: ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*##"; printf "\033[36m\033[0m"} /^[a-zA-Z0-9_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)


check_sprocket:
	@echo "Checking your sprocket version..."
	@if ! command -v sprocket >/dev/null 2>&1; then \
		echo >&2 "Error: sprocket is not installed or not in PATH. Install sprocket (https://sprocket.bio/installation.html)"; \
		exit 1; \
	else \
		sprocket_version=$$(sprocket --version | awk '{print $$2}'); \
		if printf '%s\n%s\n' "$(SPROCKET_MIN)" "$$sprocket_version" | sort -V -C; then \
			echo "sprocket $$sprocket_version is installed and compatible"; \
		else \
			echo >&2 "Error: sprocket $$sprocket_version is older than the required minimum $(SPROCKET_MIN)."; \
			echo >&2 "Please upgrade sprocket: https://sprocket.bio/installation.html"; \
			exit 1; \
		fi; \
	fi;

check_uv:
	@echo "Checking if uv is available..."
	@if ! command -v uv >/dev/null 2>&1; then \
		echo >&2 "Error: uv is not installed or not in PATH. Install uv (https://docs.astral.sh/uv/getting-started/installation/)"; \
		exit 1; \
	else \
		echo "uv version $$(uv --version | awk '{print $$2}')"; \
		uv run --python 3.13 --with miniwdl==$(MINIWDL) python -c "from importlib.metadata import version; print(f'miniwdl v{version(\"miniwdl\")}')"; \
	fi;

check_name:
	@if [ "$(NAME)" != "*" ]; then \
		if [ ! -d "modules/$(NAME)" ] && [ ! -d "pipelines/$(NAME)" ]; then \
			echo >&2 "Error: '$(NAME)' not found in modules/ or pipelines/ directory"; \
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

lint_sprocket: check_sprocket check_name ## Run sprocket lint on modules and pipelines (use NAME=foo for specific item)
	@echo "Running sprocket lint..."
	@for dir in modules/$(NAME)/ pipelines/$(NAME)/; do \
		if [ -d "$$dir" ]; then \
			echo "Linting $$dir"; \
			sprocket lint "$$dir"; \
		fi; \
	done

lint_miniwdl: check_uv check_name ## Run miniwdl lint on modules and pipelines (use NAME=foo, VERBOSE=1 for details)
	@echo "Running miniwdl lint..."
	@for file in modules/$(NAME)/*.wdl pipelines/$(NAME)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Linting $$file"; \
			if [ "$(VERBOSE)" = "1" ]; then \
				uv run --python 3.13 --with miniwdl==$(MINIWDL) miniwdl check "$$file"; \
			else \
				uv run --python 3.13 --with miniwdl==$(MINIWDL) miniwdl check "$$file" >/dev/null; \
			fi; \
		fi; \
	done

lint_womtool: check_java check_womtool check_name ## Run WOMtool validate on modules and pipelines (use NAME=foo for specific item)
	@echo "Running WOMtool validate..."
	@set -e; for file in modules/$(NAME)/*.wdl pipelines/$(NAME)/*.wdl; do \
		if [ -f "$$file" ]; then \
			echo "Validating $$file"; \
			java -jar $(WOMTOOL_JAR) validate "$$file"; \
		fi; \
	done


lint_cirro: check_name ## Validate .cirro configurations in pipelines (use NAME=foo for specific item)
	@echo "Validating Cirro configurations..."
	@if [ "$(NAME)" != "*" ]; then \
		python3 .github/scripts/validate_cirro.py $(NAME); \
	else \
		python3 .github/scripts/validate_cirro.py; \
	fi

lint: lint_sprocket lint_miniwdl lint_womtool lint_cirro ## Run all linting checks

##@ Run

run_sprocket: check_sprocket check_name ## Run sprocket on testrun WDLs (use NAME=foo, TYPE=modules|pipelines, TARGET=ci|hpc, SPROCKET_CONFIG=path, NUM_RETRIES=N)
	@echo "Running sprocket on test-run WDLs (target: $(TARGET))..."
	@failed=""; \
	dirs=""; \
	max_attempts=$$(($(NUM_RETRIES) + 1)); \
	if [ "$(TYPE)" = "modules" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs modules/$(NAME)/"; fi; \
	if [ "$(TYPE)" = "pipelines" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs pipelines/$(NAME)/"; fi; \
	for dir in $$dirs; do \
		if [ ! -d "$$dir" ]; then continue; fi; \
		wdl=""; \
		if [ "$(TARGET)" = "hpc" ] && [ -f "$$dir/testrun_hpc.wdl" ]; then \
			wdl="$$dir/testrun_hpc.wdl"; \
		elif [ -f "$$dir/testrun.wdl" ]; then \
			wdl="$$dir/testrun.wdl"; \
		fi; \
		if [ -z "$$wdl" ]; then continue; fi; \
		name=$$(basename $$dir); \
		echo "... Running $$name ($$wdl)"; \
		entrypoint=$$(grep '^workflow ' "$$wdl" | awk '{print $$2}' | tr -d '{'); \
		echo "... Using entrypoint: $$entrypoint"; \
		attempt=1; passed=0; \
		while [ $$attempt -le $$max_attempts ]; do \
			if [ $$attempt -gt 1 ]; then echo "[retry] Attempt $$attempt/$$max_attempts for $$name"; fi; \
			if sprocket run $(SPROCKET_CONFIG_FLAG) "$$wdl" --target $$entrypoint; then \
				passed=1; break; \
			fi; \
			attempt=$$((attempt + 1)); \
		done; \
		if [ $$passed -eq 0 ]; then \
			failed="$$failed $$name"; \
			echo "... FAILED: $$name (sprocket, after $$max_attempts attempts)"; \
		elif [ $$attempt -gt 1 ]; then \
			echo "[retry] $$name passed on attempt $$attempt"; \
		fi; \
	done; \
	if [ -n "$$failed" ]; then \
		echo ""; \
		echo "=== SPROCKET FAILURES ==="; \
		echo "The following items failed with sprocket:$$failed"; \
		exit 1; \
	else \
		echo "All sprocket runs passed."; \
	fi

run_miniwdl: check_uv check_name ## Run miniwdl on testrun WDLs (use NAME=foo, TYPE=modules|pipelines, TARGET=ci|hpc)
	@echo "Running miniwdl on test-run WDLs (target: $(TARGET))..."
	@failed=""; \
	dirs=""; \
	if [ "$(TYPE)" = "modules" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs modules/$(NAME)/"; fi; \
	if [ "$(TYPE)" = "pipelines" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs pipelines/$(NAME)/"; fi; \
	for dir in $$dirs; do \
		if [ ! -d "$$dir" ]; then continue; fi; \
		wdl=""; \
		if [ "$(TARGET)" = "hpc" ] && [ -f "$$dir/testrun_hpc.wdl" ]; then \
			wdl="$$dir/testrun_hpc.wdl"; \
		elif [ -f "$$dir/testrun.wdl" ]; then \
			wdl="$$dir/testrun.wdl"; \
		fi; \
		if [ -z "$$wdl" ]; then continue; fi; \
		echo "... Running $$(basename $$dir) ($$wdl)"; \
		if ! uv run --python 3.13 --with miniwdl==$(MINIWDL) miniwdl run "$$wdl"; then \
			failed="$$failed $$(basename $$dir)"; \
			echo "... FAILED: $$(basename $$dir) (miniwdl)"; \
		fi; \
	done; \
	if [ -n "$$failed" ]; then \
		echo ""; \
		echo "=== MINIWDL FAILURES ==="; \
		echo "The following items failed with miniwdl:$$failed"; \
		exit 1; \
	else \
		echo "All miniwdl runs passed."; \
	fi

run_cromwell: check_java check_cromwell check_name ## Run Cromwell on testrun WDLs (use NAME=foo, TYPE=modules|pipelines, TARGET=ci|hpc)
	@echo "Running Cromwell on test-run WDLs (target: $(TARGET))..."
	@failed=""; \
	dirs=""; \
	if [ "$(TYPE)" = "modules" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs modules/$(NAME)/"; fi; \
	if [ "$(TYPE)" = "pipelines" ] || [ "$(TYPE)" = "all" ]; then dirs="$$dirs pipelines/$(NAME)/"; fi; \
	for dir in $$dirs; do \
		if [ ! -d "$$dir" ]; then continue; fi; \
		wdl=""; \
		if [ "$(TARGET)" = "hpc" ] && [ -f "$$dir/testrun_hpc.wdl" ]; then \
			wdl="$$dir/testrun_hpc.wdl"; \
		elif [ -f "$$dir/testrun.wdl" ]; then \
			wdl="$$dir/testrun.wdl"; \
		fi; \
		if [ -z "$$wdl" ]; then continue; fi; \
		echo "... Running $$(basename $$dir) ($$wdl)"; \
		if ! java -jar $(CROMWELL_JAR) run "$$wdl"; then \
			failed="$$failed $$(basename $$dir)"; \
			echo "... FAILED: $$(basename $$dir) (Cromwell)"; \
		fi; \
	done; \
	if [ -n "$$failed" ]; then \
		echo ""; \
		echo "=== CROMWELL FAILURES ==="; \
		echo "The following items failed with Cromwell:$$failed"; \
		exit 1; \
	else \
		echo "All Cromwell runs passed."; \
	fi

run: run_sprocket run_miniwdl run_cromwell ## Run all engines on testrun.wdl files

##@ Documentation

docs-preview: check_sprocket check_uv ## Build and serve documentation preview locally
	@echo "Building documentation preview..."
	@echo "Step 1/7: Checking for uncommitted changes..."
	@if git diff --quiet HEAD -- modules/ pipelines/ && \
	   ! git ls-files --others --exclude-standard modules/ pipelines/ | grep -q .; then \
		echo "... No uncommitted changes detected"; \
	else \
		echo ""; \
		echo "WARNING: You have uncommitted changes in modules/ or pipelines/"; \
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
	@git stash push -u -m "docs-preview: temporary stash before building docs" -- modules/ pipelines/
	@echo "Step 3/7: Creating preambles..."
	@uv run --python 3.13 .github/scripts/make_preambles.py
	@echo "Step 4/7: Creating .sprocketignore..."
	@printf '%s\n' '# Excluding test run WDL'\''s from documentation builds' 'modules/**/testrun.wdl' 'modules/**/testrun_hpc.wdl' 'pipelines/**/testrun.wdl' 'pipelines/**/testrun_hpc.wdl' > .sprocketignore
	@echo "Step 5/7: Building docs with sprocket..."
	@sprocket dev doc -v --homepage docs-README.md --logo WILDSWDLNameLogo.svg .
	@echo "Step 6/7: Post-processing documentation..."
	@uv run --python 3.13 .github/scripts/postprocess_docs.py
	@echo "Step 7/7: Restoring original state..."
	@git restore README.md modules/ pipelines/ 2>/dev/null || true
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
